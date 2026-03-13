"""
Heat Exchanger Solver — Level 2 Tube-Segment Model
Nr × Nt × N_seg segment-by-segment with T_wall iteration convergence.
"""
import math
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Literal
from .properties import RefrigerantProperties, MoistAirProperties
from .geometry import FinTubeSpec, FinTubeGeo, MCHXSpec, MCHXGeo, generate_circuits
from .correlations import (
    compute_j_factor, select_correlations, recommend_correlation,
    compute_f_factor, recommend_f_correlation, get_available_f_correlations,
    FSIDE_CORRELATIONS,
    REFSIDE_EVAP_CORRELATIONS, REFSIDE_COND_CORRELATIONS,
    h_with_transition,
    h_single_gnielinski,
)


# ============================================================
# Input Data Classes
# ============================================================

@dataclass
class SimulationInput:
    """Complete simulation input."""
    # Heat exchanger type
    hx_type: Literal["FT", "MCHX"] = "FT"
    mode: Literal["evap", "cond"] = "evap"

    # Air inlet
    T_air_in: float = 308.15    # [K] (35°C)
    RH_in: float = 0.50        # [-]
    V_air: float = 2.0          # face velocity [m/s]
    P_atm: float = 101325.0    # [Pa]

    # Refrigerant
    fluid: str = "R410A"
    T_sat: float = 280.15      # [K] (7°C for evaporator)
    m_ref: float = 0.02        # total mass flow rate [kg/s]
    x_in: float = 0.2          # inlet quality (evap) or 1.0 (cond)
    T_ref_in: Optional[float] = None  # [K] actual inlet temp (for single-phase entry)
    superheat: float = 5.0     # [K] target superheat (evap)
    subcool: float = 5.0       # [K] target subcool (cond)

    # Geometry (FT)
    ft_spec: Optional[FinTubeSpec] = None
    # Geometry (MCHX)
    mchx_spec: Optional[MCHXSpec] = None

    # Flow arrangement
    flow_arrangement: Literal["counter", "parallel"] = "counter"

    # Solver — inner (T_wall per segment)
    alpha: float = 0.7         # under-relaxation
    max_iter: int = 12         # T_wall iteration limit
    tol_T: float = 0.05        # [K]
    tol_Q: float = 0.5         # [W]

    # Solver — outer (air-refrigerant coupling)
    max_outer: int = 30        # outer iteration limit
    outer_tol_pct: float = 0.1 # [%] relative Q convergence
    outer_tol_T: float = 0.1   # [K] T_air_out convergence


@dataclass
class SegmentResult:
    """Result for a single segment."""
    row: int = 0
    tube: int = 0
    seg: int = 0
    Q: float = 0.0             # heat transfer [W]
    Q_sensible: float = 0.0
    Q_latent: float = 0.0
    T_wall: float = 0.0        # [K]
    T_air_local: float = 0.0   # [K]
    T_ref: float = 0.0         # [K]
    x_ref: float = 0.0         # quality
    h_i: float = 0.0           # refrigerant-side HTC
    h_o: float = 0.0           # air-side HTC
    eta_o: float = 0.0         # overall fin efficiency
    is_wet: bool = False
    converged: bool = False
    n_iter: int = 0


@dataclass
class SimulationResult:
    """Complete simulation result."""
    Q_total: float = 0.0
    Q_sensible: float = 0.0
    Q_latent: float = 0.0
    SHR: float = 0.0
    T_air_out: float = 0.0
    W_air_out: float = 0.0
    RH_out: float = 0.0
    x_ref_out: float = 0.0
    T_ref_out: float = 0.0
    dp_air: float = 0.0
    segments: List[SegmentResult] = field(default_factory=list)
    row_Q: List[float] = field(default_factory=list)
    correlations_used: Dict = field(default_factory=dict)
    convergence: Dict = field(default_factory=dict)
    error: str = ""


# ============================================================
# Main Solver
# ============================================================

class HXSolver:
    """Level 2 Tube-Segment heat exchanger solver."""

    def __init__(self, inp: SimulationInput):
        self.inp = inp
        self.ref = RefrigerantProperties(inp.fluid)
        self.air = MoistAirProperties()

        # Build geometry
        if inp.hx_type == "FT":
            if inp.ft_spec is None:
                inp.ft_spec = FinTubeSpec()
            self.geo = FinTubeGeo.from_spec(inp.ft_spec)
            self.spec = inp.ft_spec
            self.Di = inp.ft_spec.Di
            self.Nr = inp.ft_spec.Nr
            self.Nt = inp.ft_spec.Nt
            self.N_seg = inp.ft_spec.N_seg
        else:
            if inp.mchx_spec is None:
                inp.mchx_spec = MCHXSpec()
            self.geo = MCHXGeo.from_spec(inp.mchx_spec)
            self.spec = inp.mchx_spec
            self.Di = self.geo.Dh_ref
            self.Nr = inp.mchx_spec.n_slabs
            self.Nt = inp.mchx_spec.N_tubes
            self.N_seg = inp.mchx_spec.N_seg

        # Auto-select correlations
        fin_type = inp.ft_spec.fin_type if inp.hx_type == "FT" else "louver"
        Pt = inp.ft_spec.Pt if inp.hx_type == "FT" else 0.01
        Pl = inp.ft_spec.Pl if inp.hx_type == "FT" else 0.01
        self.corr = select_correlations(inp.hx_type, self.Di, fin_type, Pt, Pl)

    def solve(self) -> SimulationResult:
        """Run the full simulation."""
        inp = self.inp
        result = SimulationResult()
        result.correlations_used = self.corr

        try:
            return self._solve_internal()
        except Exception as e:
            result.error = str(e)
            return result

    def _solve_internal(self) -> SimulationResult:
        inp = self.inp
        ref = self.ref
        P_sat = ref.P_sat(inp.T_sat)

        # Air inlet state
        W_in = self.air.W_from_TRH(inp.T_air_in, inp.RH_in, inp.P_atm)
        h_air_in = self.air.h_simple(inp.T_air_in, W_in)
        T_dp = self.air.Tdp_from_TW(inp.T_air_in, W_in, inp.P_atm)

        # Air mass flow rate
        rho_air = self.air.rho_air(inp.T_air_in, W_in, inp.P_atm)
        A_fr = self.geo.A_fr
        m_air = rho_air * inp.V_air * A_fr

        # ── Build circuit paths ──
        if inp.hx_type == "FT":
            spec = inp.ft_spec
            if spec.circuits and spec.circuit_mode == "custom":
                circuits = spec.circuits
            else:
                circuits = generate_circuits(
                    self.Nr, self.Nt, spec.circuit_mode, inp.flow_arrangement)
        else:
            # MCHX: all parallel (each tube is one circuit through all slabs)
            circuits = generate_circuits(
                self.Nr, self.Nt, "row_parallel", inp.flow_arrangement)

        n_circ = len(circuits)

        # Refrigerant per circuit
        m_ref_circ = inp.m_ref / n_circ
        if inp.hx_type == "FT":
            A_cs_ref = math.pi * self.Di ** 2 / 4
            G_ref = m_ref_circ / A_cs_ref if A_cs_ref > 0 else 100
        else:
            A_cs_port = inp.mchx_spec.ch_width * inp.mchx_spec.ch_height
            m_per_port = m_ref_circ / inp.mchx_spec.n_ports
            G_ref = m_per_port / A_cs_port if A_cs_port > 0 else 100

        # Air mass flux & h_o
        G_air = m_air / self.geo.A_c if self.geo.A_c > 0 else 5.0
        h_o = self._compute_h_o(G_air, inp.T_air_in)
        h_fg = ref.h_fg(P_sat)

        # ── Outer iteration for air temperature convergence ──
        # Air state per column per row: T_air_2d[col][row]
        # Each column has INDEPENDENT air stream flowing through rows 0→Nr-1
        Ns = self.N_seg
        T_air_3d = [[[inp.T_air_in] * (self.Nr + 1) for _ in range(Ns)]
                     for _ in range(self.Nt)]
        W_air_3d = [[[W_in] * (self.Nr + 1) for _ in range(Ns)]
                     for _ in range(self.Nt)]
        m_air_cell = m_air / max(self.Nt * Ns, 1)

        seg_dict = {}
        max_outer = inp.max_outer
        outer_tol_pct = inp.outer_tol_pct
        outer_tol_T = inp.outer_tol_T
        Q_prev_outer = 0.0
        T_air_out_prev = 0.0
        outer_converged = False
        outer_history = []
        omega = 1.0  # under-relaxation: adaptive

        for outer_iter in range(max_outer):
            seg_dict.clear()
            circ_outlets = []

            # ── Process each circuit ──
            for circ_idx, path in enumerate(circuits):
                x_ref = inp.x_in
                T_ref = inp.T_sat
                if inp.T_ref_in is not None and (inp.x_in >= 1.0 or inp.x_in <= 0.0):
                    T_ref = inp.T_ref_in

                for pass_idx, (row_idx, col_idx) in enumerate(path):
                    # Alternate segment direction per tube pass (U-bend)
                    # Even pass: seg 0 → N-1, Odd pass: seg N-1 → 0
                    if pass_idx % 2 == 0:
                        seg_order = range(Ns)
                    else:
                        seg_order = range(Ns - 1, -1, -1)

                    for seg_idx in seg_order:
                        T_air_local = T_air_3d[col_idx][seg_idx][row_idx]
                        W_air_local = W_air_3d[col_idx][seg_idx][row_idx]

                        seg_result = self._solve_segment(
                            row=row_idx, tube=col_idx, seg=seg_idx,
                            T_air=T_air_local, W_air=W_air_local, T_dp=T_dp,
                            x_ref=x_ref, T_ref=T_ref,
                            P_sat=P_sat, G_ref=G_ref, G_air=G_air,
                            h_o=h_o, m_ref_tube=m_ref_circ,
                        )
                        seg_dict[(col_idx, row_idx, seg_idx)] = seg_result

                        # Update refrigerant state
                        if 0 <= x_ref <= 1:
                            if h_fg > 0 and m_ref_circ > 0:
                                dx = seg_result.Q / (m_ref_circ * h_fg)
                                if inp.mode == "evap":
                                    x_ref += dx
                                else:
                                    x_ref -= dx
                        else:
                            if x_ref > 1:
                                try:
                                    cp_v = ref.cp_v(P_sat)
                                    if m_ref_circ > 0:
                                        T_ref += seg_result.Q / (m_ref_circ * cp_v)
                                except: pass
                            elif x_ref < 0:
                                try:
                                    cp_l = ref.cp_l(P_sat)
                                    if m_ref_circ > 0:
                                        T_ref -= seg_result.Q / (m_ref_circ * cp_l)
                                except: pass

                circ_outlets.append((x_ref, T_ref))

            # ── Check outer convergence ──
            Q_this = sum(s.Q for s in seg_dict.values())
            T_air_out_this = sum(T_air_3d[c][s][self.Nr]
                                 for c in range(self.Nt) for s in range(Ns)) / max(self.Nt * Ns, 1)
            dQ_outer = abs(Q_this - Q_prev_outer)
            dT_outer = abs(T_air_out_this - T_air_out_prev)
            rel_dQ = dQ_outer / max(abs(Q_this), 1.0) * 100

            outer_history.append({
                "iter": outer_iter + 1,
                "Q": round(Q_this, 2),
                "dQ": round(dQ_outer, 2),
                "dQ_pct": round(rel_dQ, 3),
                "T_air_out": round(T_air_out_this - 273.15, 2),
                "dT": round(dT_outer, 3),
                "omega": round(omega, 3),
            })

            if outer_iter > 0 and rel_dQ < outer_tol_pct and dT_outer < outer_tol_T:
                outer_converged = True
                Q_prev_outer = Q_this
                T_air_out_prev = T_air_out_this
                break

            # ── Adaptive under-relaxation ──
            # Detect oscillation: if Q flipped sign between last 2 iterations
            if outer_iter >= 2:
                h1 = outer_history[-2]
                h2 = outer_history[-3]
                dQ_prev = h1["Q"] - h2["Q"]
                dQ_curr = Q_this - h1["Q"]
                if dQ_prev * dQ_curr < 0:
                    # Oscillation detected → reduce omega
                    omega = max(omega * 0.7, 0.3)
                else:
                    # Monotone → increase omega toward 1.0
                    omega = min(omega * 1.1, 0.9)
            elif outer_iter == 1:
                omega = 0.6  # conservative start after first full update

            Q_prev_outer = Q_this
            T_air_out_prev = T_air_out_this

            # ── Update per-(col, seg) air pencil with adaptive omega ──
            h_fg_water = 2501000.0
            for col_idx in range(self.Nt):
                for seg_idx in range(Ns):
                    T_air_3d[col_idx][seg_idx][0] = inp.T_air_in
                    W_air_3d[col_idx][seg_idx][0] = W_in
                    for row_idx in range(self.Nr):
                        sr = seg_dict.get((col_idx, row_idx, seg_idx))
                        Q_seg = sr.Q if sr else 0.0
                        Q_lat = sr.Q_latent if sr else 0.0
                        if m_air_cell > 0:
                            h_cur = self.air.h_simple(
                                T_air_3d[col_idx][seg_idx][row_idx],
                                W_air_3d[col_idx][seg_idx][row_idx])
                            if inp.mode == "evap":
                                h_next = h_cur - Q_seg / m_air_cell
                            else:
                                h_next = h_cur + Q_seg / m_air_cell
                            W_rem = Q_lat / (m_air_cell * h_fg_water)
                            W_next = max(W_air_3d[col_idx][seg_idx][row_idx] - W_rem, 0)
                            T_calc = self.air.T_from_h_simple(h_next, W_next)
                            # Under-relax: blend new with old
                            T_old = T_air_3d[col_idx][seg_idx][row_idx + 1]
                            T_air_3d[col_idx][seg_idx][row_idx + 1] = omega * T_calc + (1 - omega) * T_old
                            W_air_3d[col_idx][seg_idx][row_idx + 1] = omega * W_next + (1 - omega) * W_air_3d[col_idx][seg_idx][row_idx + 1]
                        else:
                            T_air_3d[col_idx][seg_idx][row_idx + 1] = T_air_3d[col_idx][seg_idx][row_idx]
                            W_air_3d[col_idx][seg_idx][row_idx + 1] = W_air_3d[col_idx][seg_idx][row_idx]

        # ── Compile results ──
        all_segments = []
        for col_idx in range(self.Nt):
            for row_idx in range(self.Nr):
                for seg_idx in range(Ns):
                    key = (col_idx, row_idx, seg_idx)
                    if key in seg_dict:
                        all_segments.append(seg_dict[key])

        Q_total = sum(s.Q for s in all_segments)
        Q_lat_total = sum(s.Q_latent for s in all_segments)
        Q_sen = Q_total - Q_lat_total

        x_ref_out_avg = sum(o[0] for o in circ_outlets) / n_circ if circ_outlets else inp.x_in
        T_ref_out_avg = sum(o[1] for o in circ_outlets) / n_circ if circ_outlets else inp.T_sat

        T_air_out = sum(T_air_3d[c][s][self.Nr]
                        for c in range(self.Nt) for s in range(Ns)) / max(self.Nt * Ns, 1)
        W_air_out = sum(W_air_3d[c][s][self.Nr]
                        for c in range(self.Nt) for s in range(Ns)) / max(self.Nt * Ns, 1)

        row_Q = [sum(seg_dict.get((t, r, s), SegmentResult()).Q
                     for t in range(self.Nt) for s in range(Ns))
                 for r in range(self.Nr)]

        result = SimulationResult()
        result.Q_total = Q_total
        result.Q_sensible = Q_sen
        result.Q_latent = Q_lat_total
        result.SHR = Q_sen / Q_total if Q_total > 0 else 1.0
        result.T_air_out = T_air_out
        result.W_air_out = W_air_out
        result.x_ref_out = x_ref_out_avg
        result.T_ref_out = T_ref_out_avg
        result.segments = all_segments
        result.row_Q = row_Q
        result.correlations_used = self.corr
        result.convergence = {
            "outer_converged": outer_converged,
            "outer_iterations": outer_iter + 1,
            "outer_max": max_outer,
            "outer_tol_pct": outer_tol_pct,
            "outer_tol_T": outer_tol_T,
            "final_dQ": outer_history[-1]["dQ"] if outer_history else 0,
            "final_dQ_pct": outer_history[-1]["dQ_pct"] if outer_history else 0,
            "final_dT": outer_history[-1]["dT"] if outer_history else 0,
            "history": outer_history,
            "seg_converged_pct": round(
                sum(1 for s in all_segments if s.converged) / max(len(all_segments), 1) * 100, 1),
        }

        # Store circuit info
        self.corr["n_circuits"] = n_circ
        self.corr["circuit_mode"] = inp.ft_spec.circuit_mode if inp.hx_type == "FT" else "row_parallel"

        # Air-side pressure drop
        result.dp_air = self._compute_dp_air(G_air, inp.T_air_in, T_air_out)

        # Outlet RH
        try:
            import CoolProp.CoolProp as CP
            result.RH_out = CP.HAPropsSI("R", "T", T_air_out, "W", W_air_out, "P", inp.P_atm)
        except:
            result.RH_out = 0.0

        return result

    def _solve_segment(self, row, tube, seg, T_air, W_air, T_dp,
                       x_ref, T_ref, P_sat, G_ref, G_air, h_o,
                       m_ref_tube) -> SegmentResult:
        """Solve single segment with T_wall iteration."""
        inp = self.inp
        ref = self.ref
        sr = SegmentResult(row=row, tube=tube, seg=seg)

        # Segment areas
        if inp.hx_type == "FT":
            A_i_seg = self.geo.A_i_seg
            A_o_seg = self.geo.A_total / (self.Nr * self.Nt * self.N_seg)
        else:
            A_i_seg = self.geo.A_i_seg
            A_o_seg = self.geo.A_total / (self.Nr * self.Nt * self.N_seg)

        # Initial T_wall guess
        T_w = (T_air + (T_ref if 0 <= x_ref <= 1 else T_ref)) / 2

        # Check wet/dry
        is_wet = (T_w < T_dp) and (inp.mode == "evap")

        alpha = inp.alpha
        Q_prev = 0

        for iteration in range(inp.max_iter):
            # --- Refrigerant side h_i ---
            # q_flux estimate for Kim&Mudawar
            q_flux_est = abs(T_air - T_w) * h_o if h_o > 0 else 5000

            h_i = h_with_transition(
                x=x_ref, G=G_ref, Di=self.Di, q_flux=q_flux_est,
                ref=ref, P=P_sat, mode=inp.mode,
                hx_type=inp.hx_type,
                evap_corr=self.corr.get("evap"),
                cond_corr=self.corr.get("cond"),
            )

            # --- Fin efficiency ---
            if is_wet:
                b = self._compute_b_factor(T_w, T_air, W_air, h_o)
                if inp.hx_type == "FT":
                    _, eta_o = self.geo.fin_efficiency_wet(h_o, b)
                else:
                    _, eta_o = self.geo.fin_efficiency_wet_straight(h_o, b)
            else:
                if inp.hx_type == "FT":
                    _, eta_o = self.geo.fin_efficiency_schmidt(h_o)
                else:
                    _, eta_o = self.geo.fin_efficiency_straight(h_o)

            # --- Thermal resistances ---
            R_o = 1.0 / (eta_o * h_o * A_o_seg) if (eta_o * h_o * A_o_seg) > 0 else 1e6
            R_i = 1.0 / (h_i * A_i_seg) if (h_i * A_i_seg) > 0 else 1e6
            UA = 1.0 / (R_o + R_i)

            # --- NTU-epsilon ---
            C_air_seg = 1006.0 * (m_ref_tube * self.Nt)  # approximate
            # For segment: use simple Q = UA × LMTD approach
            if inp.mode == "evap":
                T_ref_eff = T_ref if x_ref > 1 else inp.T_sat
                dT = T_air - T_ref_eff
            else:
                T_ref_eff = T_ref if x_ref < 0 else inp.T_sat
                dT = T_ref_eff - T_air

            Q_air = UA * dT if dT > 0 else 0

            # --- Wet surface: additional latent heat ---
            Q_lat = 0
            if is_wet and inp.mode == "evap":
                Ws_wall = self.air.Ws_from_T(T_w, inp.P_atm)
                if W_air > Ws_wall:
                    h_fg_w = 2501000.0
                    Q_lat = eta_o * h_o * A_o_seg * h_fg_w * (W_air - Ws_wall) / 1006.0
                    Q_lat = max(Q_lat, 0)

            Q_total_seg = Q_air + Q_lat

            # --- Update T_wall ---
            T_w_new = T_ref_eff + Q_total_seg * R_i if inp.mode == "evap" else T_ref_eff - Q_total_seg * R_i
            T_w_calc = alpha * T_w_new + (1 - alpha) * T_w

            # Convergence check
            dT_w = abs(T_w_calc - T_w)
            dQ = abs(Q_total_seg - Q_prev)

            T_w = T_w_calc
            Q_prev = Q_total_seg

            # Re-check wet condition
            is_wet = (T_w < T_dp) and (inp.mode == "evap")

            if dT_w < inp.tol_T and dQ < inp.tol_Q:
                sr.converged = True
                sr.n_iter = iteration + 1
                break

        sr.Q = Q_total_seg
        sr.Q_sensible = Q_air
        sr.Q_latent = Q_lat
        sr.T_wall = T_w
        sr.T_air_local = T_air
        sr.T_ref = T_ref_eff
        sr.x_ref = x_ref
        sr.h_i = h_i
        sr.h_o = h_o
        sr.eta_o = eta_o
        sr.is_wet = is_wet

        return sr

    def _compute_h_o(self, G_air: float, T_air: float) -> float:
        """Compute air-side HTC using selected j-factor correlation."""
        inp = self.inp
        mu = self.air.mu_air(T_air)
        Pr = self.air.Pr_air(T_air)
        cp = 1006.0

        corr_id = self.corr["air_j"]

        if inp.hx_type == "FT":
            spec = inp.ft_spec
            Dc = self.geo.Dc
            Re_Dc = G_air * Dc / mu

            j = compute_j_factor(
                corr_id,
                Re_Dc=Re_Dc, Nr=spec.Nr, Dc=Dc,
                Pt=spec.Pt, Pl=spec.Pl, FPI=spec.FPI,
                fin_thickness=spec.fin_thickness,
                # Wavy params
                Xa=spec.wavy_amplitude, wave_length=spec.wavy_wavelength,
                # Louver params
                Lp=spec.louver_pitch, theta=spec.louver_angle,
                # Slit params
                slit_height=spec.slit_height, slit_width=spec.slit_width,
                n_slits=spec.n_slits,
            )
        else:  # MCHX
            spec = inp.mchx_spec
            Re_Lp = G_air * spec.louver_pitch / mu
            j = compute_j_factor(
                corr_id,
                Re_Lp=Re_Lp, Lp=spec.louver_pitch,
                theta=spec.louver_angle, Fp=spec.fin_pitch,
                fin_thickness=spec.fin_thickness,
                Fl=spec.fin_height, Td=spec.D,
            )

        h_o = j * G_air * cp / Pr ** (2 / 3)
        return max(h_o, 10.0)

    def _compute_b_factor(self, T_wall: float, T_air: float, W_air: float,
                          h_o: float) -> float:
        """
        Compute b-factor for wet surface at T_fin_avg.
        b = cp_s/cp_a = 1 + h_fg × (dWs/dT) / cp_air
        Iterative: T_fin_avg → b → η_wet → T_fin_avg (3 iterations)
        """
        cp_air = 1006.0
        h_fg = 2501000.0

        # Initial dry fin efficiency for T_fin estimate
        if self.inp.hx_type == "FT":
            _, eta_dry = self.geo.fin_efficiency_schmidt(h_o)
        else:
            _, eta_dry = self.geo.fin_efficiency_straight(h_o)

        T_fin = T_air - eta_dry * (T_air - T_wall)

        for _ in range(3):
            T_fin = max(T_fin, 250.0)
            T_fin = min(T_fin, 370.0)
            dWs_dT = self.air.dWs_dT(T_fin)
            b = 1.0 + h_fg * dWs_dT / cp_air
            b = max(b, 1.0)

            # Recompute fin temp with wet efficiency
            if self.inp.hx_type == "FT":
                _, eta_wet = self.geo.fin_efficiency_wet(h_o, b)
            else:
                _, eta_wet = self.geo.fin_efficiency_wet_straight(h_o, b)

            T_fin = T_air - eta_wet * (T_air - T_wall)

        return b

    def _compute_dp_air(self, G_air: float, T_in: float, T_out: float) -> float:
        """
        Air-side pressure drop — Kays & London (1984) formulation.
        Uses f-factor from registry (auto or manual selection).
        """
        inp = self.inp
        T_avg = (T_in + T_out) / 2
        mu = self.air.mu_air(T_avg)
        rho_in = self.air.rho_air(T_in, 0.01)
        rho_out = self.air.rho_air(T_out, 0.01)
        rho_m = (rho_in + rho_out) / 2

        if inp.hx_type == "FT":
            spec = inp.ft_spec
            Dc = self.geo.Dc
            Re_Dc = G_air * Dc / mu
            sigma = self.geo.sigma

            # Determine f-factor correlation
            f_corr_id = getattr(self, '_f_corr_id', None)
            if not f_corr_id:
                # Auto-select default based on fin type
                fin_type_map = {
                    "plain": "f_wang2000_plain",
                    "wavy": "f_wang1999_wavy",
                    "slit": "f_wang2001_slit",
                    "louver": "f_louver_enhanced",
                }
                f_corr_id = fin_type_map.get(spec.fin_type, "f_wang2000_plain")

            # Build kwargs for the correlation
            f_kwargs = dict(
                Re_Dc=Re_Dc, Nr=spec.Nr, Dc=Dc,
                Pt=spec.Pt, Pl=spec.Pl, FPI=spec.FPI,
                fin_thickness=spec.fin_thickness,
                Xa=spec.wavy_amplitude, wave_length=spec.wavy_wavelength,
                slit_height=spec.slit_height, slit_width=spec.slit_width,
                n_slits=spec.n_slits,
                Lp=spec.louver_pitch, theta=spec.louver_angle,
                Ao_Ac=self.geo.A_total / self.geo.A_c if self.geo.A_c > 0 else 200,
            )

            try:
                f = compute_f_factor(f_corr_id, **f_kwargs)
            except Exception:
                # Fallback to plain
                f = compute_f_factor("f_wang2000_plain", **f_kwargs)

            # Store used f-correlation info
            self.corr["air_f"] = f_corr_id
            self.corr["f_value"] = round(f, 6)
            self.corr["Re_Dc_f"] = round(Re_Dc, 1)

            A_ratio = self.geo.A_total / self.geo.A_c if self.geo.A_c > 0 else 10
            Kc = 0.42 * (1 - sigma ** 2)
            Ke = max((1 - sigma ** 2 - 0.4 * (1 - sigma ** 2) ** 1.25), 0.0)

        else:
            spec = inp.mchx_spec
            Re_Lp = G_air * spec.louver_pitch / mu

            f_corr_id = getattr(self, '_f_corr_id', None)
            if not f_corr_id:
                f_corr_id = "f_chang_wang1997_mchx"

            f_kwargs = dict(
                Re_Lp=Re_Lp, Lp=spec.louver_pitch,
                theta=spec.louver_angle, Fp=spec.fin_pitch,
            )

            try:
                f = compute_f_factor(f_corr_id, **f_kwargs)
            except Exception:
                f = compute_f_factor("f_chang_wang1997_mchx", **f_kwargs)

            self.corr["air_f"] = f_corr_id
            self.corr["f_value"] = round(f, 6)

            sigma = self.geo.sigma
            A_ratio = self.geo.A_total / self.geo.A_c if self.geo.A_c > 0 else 10
            Kc = 0.42 * (1 - sigma ** 2)
            Ke = max((1 - sigma ** 2 - 0.4 * (1 - sigma ** 2) ** 1.25), 0.0)

        # Kays & London full pressure drop
        dp = G_air ** 2 / (2 * rho_in) * (
            Kc +
            (1 + sigma ** 2) * (rho_in / rho_out - 1) +
            f * A_ratio * (rho_in / rho_m) -
            Ke * (rho_in / rho_out)
        )
        return max(dp, 0.0)
