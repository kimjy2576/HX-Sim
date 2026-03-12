"""
Heat Exchanger Solver — Level 2 Tube-Segment Model
Nr × Nt × N_seg segment-by-segment with T_wall iteration convergence.
"""
import math
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Literal
from .properties import RefrigerantProperties, MoistAirProperties
from .geometry import FinTubeSpec, FinTubeGeo, MCHXSpec, MCHXGeo
from .correlations import (
    compute_j_factor, select_correlations, recommend_correlation,
    compute_f_factor, recommend_f_correlation, get_available_f_correlations,
    FSIDE_CORRELATIONS,
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
    superheat: float = 5.0     # [K] target superheat (evap)
    subcool: float = 5.0       # [K] target subcool (cond)

    # Geometry (FT)
    ft_spec: Optional[FinTubeSpec] = None
    # Geometry (MCHX)
    mchx_spec: Optional[MCHXSpec] = None

    # Flow arrangement
    flow_arrangement: Literal["counter", "parallel"] = "counter"

    # Solver
    alpha: float = 0.7         # under-relaxation
    max_iter: int = 12         # T_wall iteration limit
    tol_T: float = 0.05        # [K]
    tol_Q: float = 0.5         # [W]


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
        if inp.hx_type == "FT":
            A_fr = self.geo.A_fr
            m_air = rho_air * inp.V_air * A_fr
        else:
            A_fr = self.geo.A_fr
            m_air = rho_air * inp.V_air * A_fr

        # Refrigerant per tube (for state tracking: dx = Q/(m_tube × h_fg))
        if inp.hx_type == "FT":
            m_ref_tube = inp.m_ref / self.Nt
        else:
            # MCHX: each tube has n_ports parallel channels sharing the tube flow
            m_ref_tube = inp.m_ref / self.Nt  # per flat tube

        # Refrigerant mass flux (for correlation calculations)
        if inp.hx_type == "FT":
            A_cs_ref = math.pi * self.Di ** 2 / 4
            G_ref = m_ref_tube / A_cs_ref if A_cs_ref > 0 else 100
        else:
            # MCHX: mass flux per individual port channel
            A_cs_port = inp.mchx_spec.ch_width * inp.mchx_spec.ch_height
            m_per_port = m_ref_tube / inp.mchx_spec.n_ports
            G_ref = m_per_port / A_cs_port if A_cs_port > 0 else 100

        # Air mass flux
        if inp.hx_type == "FT":
            G_air = m_air / self.geo.A_c if self.geo.A_c > 0 else 5.0
        else:
            G_air = m_air / self.geo.A_c if self.geo.A_c > 0 else 5.0

        # Compute air-side h_o (constant for all segments)
        h_o = self._compute_h_o(G_air, inp.T_air_in)

        # Initialize segment results
        all_segments = []
        row_Q = []

        # Air state tracking per row
        T_air_row = inp.T_air_in
        W_air_row = W_in
        h_air_row = h_air_in

        # Per-tube refrigerant state tracking
        # Each tube has independent refrigerant path through all Nr rows
        x_ref_tube = [inp.x_in] * self.Nt
        T_ref_tube = [inp.T_sat] * self.Nt

        h_fg = ref.h_fg(P_sat)

        # Row order: air always goes 0 → Nr-1
        # Counter: refrigerant enters at air row (Nr-1), exits at air row 0
        # Parallel: refrigerant enters at air row 0
        # We iterate air rows in order; for counter flow, the ref state at each
        # air row depends on what happened in later air rows. We use a simple
        # forward pass (row-by-row in air direction) — this is an approximation
        # but standard for row-by-row models.

        # ===== ROW LOOP (air direction: Row 0 → Row Nr-1) =====
        for air_row_idx in range(self.Nr):
            Q_row = 0.0
            Q_row_lat = 0.0

            # Air state for this row
            T_air_local = T_air_row
            W_air_local = W_air_row

            # ===== TUBE LOOP =====
            for tube_idx in range(self.Nt):
                x_ref = x_ref_tube[tube_idx]
                T_ref = T_ref_tube[tube_idx]

                # ===== SEGMENT LOOP =====
                for seg_idx in range(self.N_seg):
                    seg_result = self._solve_segment(
                        row=air_row_idx,
                        tube=tube_idx,
                        seg=seg_idx,
                        T_air=T_air_local,
                        W_air=W_air_local,
                        T_dp=T_dp,
                        x_ref=x_ref,
                        T_ref=T_ref,
                        P_sat=P_sat,
                        G_ref=G_ref,
                        G_air=G_air,
                        h_o=h_o,
                        m_ref_tube=m_ref_tube,
                    )
                    all_segments.append(seg_result)
                    Q_row += seg_result.Q
                    Q_row_lat += seg_result.Q_latent

                    # Update per-tube refrigerant state
                    if 0 <= x_ref <= 1:
                        if h_fg > 0 and m_ref_tube > 0:
                            if inp.mode == "evap":
                                dx = seg_result.Q / (m_ref_tube * h_fg)
                                x_ref += dx
                            else:
                                dx = seg_result.Q / (m_ref_tube * h_fg)
                                x_ref -= dx
                    else:
                        if x_ref > 1:
                            try:
                                cp_v = ref.cp_v(P_sat)
                                if m_ref_tube > 0:
                                    T_ref += seg_result.Q / (m_ref_tube * cp_v)
                            except:
                                pass
                        elif x_ref < 0:
                            try:
                                cp_l = ref.cp_l(P_sat)
                                if m_ref_tube > 0:
                                    T_ref -= seg_result.Q / (m_ref_tube * cp_l)
                            except:
                                pass

                # Store updated state for this tube
                x_ref_tube[tube_idx] = x_ref
                T_ref_tube[tube_idx] = T_ref

            # Update air state after this row (enthalpy-based)
            if m_air > 0:
                if inp.mode == "evap":
                    # Evaporator: air loses heat
                    h_air_out = h_air_row - Q_row / m_air
                else:
                    # Condenser: air gains heat
                    h_air_out = h_air_row + Q_row / m_air

                # Update W (humidity ratio) — only changes for wet evaporator
                h_fg_water = 2501000.0
                W_removed = Q_row_lat / (m_air * h_fg_water) if m_air > 0 else 0
                W_air_next = W_air_row - W_removed
                W_air_next = max(W_air_next, 0)

                # Update T from enthalpy
                T_air_next = self.air.T_from_h_simple(h_air_out, W_air_next)

                T_air_row = T_air_next
                W_air_row = W_air_next
                h_air_row = h_air_out

            row_Q.append(Q_row)

        # ===== Compile results =====
        Q_total = sum(s.Q for s in all_segments)
        Q_lat = sum(s.Q_latent for s in all_segments)
        Q_sen = Q_total - Q_lat

        result = SimulationResult()
        result.Q_total = Q_total
        result.Q_sensible = Q_sen
        result.Q_latent = Q_lat
        result.SHR = Q_sen / Q_total if Q_total > 0 else 1.0
        result.T_air_out = T_air_row
        result.W_air_out = W_air_row
        result.x_ref_out = sum(x_ref_tube) / len(x_ref_tube)
        result.T_ref_out = sum(T_ref_tube) / len(T_ref_tube)
        result.segments = all_segments
        result.row_Q = row_Q
        result.correlations_used = self.corr

        # Air-side pressure drop
        result.dp_air = self._compute_dp_air(G_air, inp.T_air_in, T_air_row)

        # Outlet RH
        try:
            import CoolProp.CoolProp as CP
            result.RH_out = CP.HAPropsSI("R", "T", T_air_row, "W", W_air_row, "P", inp.P_atm)
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
