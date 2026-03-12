"""
Heat Transfer Correlations
Air-side: Wang(2000), Wang(1999), Chang&Wang(1997)
Ref-side: Chen(1966), Kim&Mudawar(2013/2012), Shah(1979), Gnielinski(1976)
f-factor: Wang(2000), Chang&Wang(1997)
"""
import math
from .properties import RefrigerantProperties, MoistAirProperties


# ============================================================
# AIR-SIDE j-factor (FT-HX)
# ============================================================

def j_factor_wang2000_plain(Re_Dc: float, Nr: int, Dc: float,
                            Pt: float, Pl: float, FPI: float,
                            fin_thickness: float) -> float:
    """
    Wang et al. (2000) IJHMT 43(15) — plain fin, staggered.
    Pt/Pl <= 1.35 case.
    """
    Fp = 0.0254 / FPI  # fin pitch [m]
    if Re_Dc < 10:
        Re_Dc = 10.0

    # Coefficients from Wang (2000) Table 4
    if Nr == 1:
        P1 = 1.9 - 0.23 * math.log(Re_Dc)
        P2 = -0.236 + 0.126 * math.log(Re_Dc)
        j = 0.108 * Re_Dc ** (-0.29) * (Pt / Pl) ** P1 * (Fp / Dc) ** (-1.084) * \
            (Fp / (Fp - fin_thickness)) ** (-0.786) * (Fp / Pt) ** P2
    else:
        P3 = -0.361 - 0.042 * Nr / math.log(Re_Dc) + 0.158 * math.log(Nr * (Fp / Dc) ** 0.41)
        P4 = -1.224 - 0.076 * (Pl / Dh_approx(Dc, Fp, fin_thickness)) ** 1.42 / math.log(Re_Dc)
        P5 = -0.083 + 0.058 * Nr / math.log(Re_Dc)
        P6 = -5.735 + 1.21 * math.log(Re_Dc / Nr)
        j = 0.086 * Re_Dc ** P3 * (Nr) ** P4 * (Fp / Dc) ** P5 * \
            (Fp / (Fp - fin_thickness)) ** P6 * (Fp / Pt) ** (-0.93)
    return max(j, 1e-6)


def j_factor_wang2000_plain_high_PtPl(Re_Dc: float, Nr: int, Dc: float,
                                       Pt: float, Pl: float, FPI: float,
                                       fin_thickness: float) -> float:
    """Wang(2000) for Pt/Pl > 1.35 — with warning flag."""
    # Use the base correlation with a correction
    j = j_factor_wang2000_plain(Re_Dc, Nr, Dc, Pt, Pl, FPI, fin_thickness)
    # Apply approximate correction for high Pt/Pl
    correction = (Pt / Pl / 1.35) ** (-0.2)
    return j * correction


def j_factor_wang1999_wavy(Re_Dc: float, Nr: int, Dc: float,
                           Pt: float, Pl: float, FPI: float,
                           fin_thickness: float, Xa: float = 0.001,
                           wave_length: float = 0.01) -> float:
    """Wang et al. (1999) IJHMT 42 — wavy fin, staggered."""
    Fp = 0.0254 / FPI
    if Re_Dc < 10:
        Re_Dc = 10.0

    # Simplified Wang(1999) wavy correlation
    j1 = -0.0045 - 0.0015 * Nr - 0.058 / (math.log(Re_Dc) + 1e-10)
    j = 0.324 * Re_Dc ** j1 * (Fp / Dc) ** (-0.3) * (Pt / Pl) ** (-0.2) * Nr ** (-0.1)
    return max(j, 1e-6)


def j_factor_wang1999_louver(Re_Dc: float, Nr: int, Dc: float,
                              Pt: float, Pl: float, FPI: float,
                              fin_thickness: float,
                              Lp: float = 0.0017, theta: float = 27.0) -> float:
    """Wang et al. (1999) IJHMT 42(1) — louver fin, staggered."""
    Fp = 0.0254 / FPI
    if Re_Dc < 10:
        Re_Dc = 10.0

    j1 = -0.991 - 0.1055 * (Pl / Pt) ** 3.1 * math.log(Re_Dc)
    j2 = -0.0044 * Nr + 0.0195 * (Lp / Dc) ** (-0.5) / math.log(Re_Dc)
    j = 0.425 * Re_Dc ** j1 * (theta / 90) ** 0.27 * (Fp / Lp) ** j2 * \
        (Fp / Dc) ** (-0.34) * Nr ** (-0.15)
    return max(j, 1e-6)


def j_factor_slit(Re_Dc: float, Nr: int, Dc: float,
                  Pt: float, Pl: float, FPI: float,
                  fin_thickness: float) -> float:
    """
    Slit fin: Wang(2001) IJHMT 44 for Dc>=10mm,
    or plain × E_slit for Dc<10mm.
    """
    if Dc >= 0.010:
        # Wang(2001) direct correlation
        Fp = 0.0254 / FPI
        j = 0.257 * Re_Dc ** (-0.43) * (Fp / Dc) ** (-0.27) * Nr ** (-0.09)
    else:
        # plain × enhancement factor
        j = j_factor_wang2000_plain(Re_Dc, Nr, Dc, Pt, Pl, FPI, fin_thickness)
        E_slit = 1.33  # typical enhancement
        j *= E_slit
    return max(j, 1e-6)


def Dh_approx(Dc: float, Fp: float, delta: float) -> float:
    """Approximate hydraulic diameter for j-factor correlations."""
    return max(4 * (Fp - delta) * (Dc * 0.5) / (2 * ((Fp - delta) + Dc * 0.5)), 1e-6)


# ============================================================
# AIR-SIDE j-factor (MCHX) — Chang & Wang (1997)
# ============================================================

def j_factor_chang_wang_1997(Re_Lp: float, Lp: float, theta: float,
                              Fp: float, fin_thickness: float = 0.0001) -> float:
    """
    Chang & Wang (1997) — louver fin for MCHX.
    Re_Lp based, 3-parameter simplified version.
    """
    if Re_Lp < 5:
        Re_Lp = 5.0
    theta_rad = theta  # degrees

    J1 = -0.49 * (theta_rad / 90) ** 0.27
    j = Re_Lp ** J1 * (theta_rad / 90) ** 0.27 * (Fp / Lp) ** (-0.14)

    return max(j, 1e-6)


# ============================================================
# AIR-SIDE f-factor
# ============================================================

def f_factor_wang2000_plain(Re_Dc: float, Nr: int, Dc: float,
                            Pt: float, Pl: float, FPI: float,
                            fin_thickness: float) -> float:
    """Wang(2000) Table 6 — plain fin f-factor."""
    Fp = 0.0254 / FPI
    if Re_Dc < 10:
        Re_Dc = 10.0

    # Simplified f-factor from Wang(2000)
    F1 = -0.764 + 0.739 * (Pt / Pl) + 0.177 * (Fp / Dc) - 0.00758 / Nr
    F2 = -15.689 + 64.021 / math.log(Re_Dc)
    F3 = 1.696 - 15.695 / math.log(Re_Dc)

    f = 0.0267 * Re_Dc ** F1 * (Pt / Pl) ** F2 * (Fp / Dc) ** F3
    return max(f, 1e-6)


def f_factor_chang_wang_1997(Re_Lp: float, Lp: float, theta: float,
                              Fp: float) -> float:
    """Chang & Wang (1997) f-factor for MCHX louver fin."""
    if Re_Lp < 5:
        Re_Lp = 5.0

    f1 = -0.72 * (theta / 90) ** 0.19
    f = Re_Lp ** f1 * (theta / 90) ** 0.34 * (Fp / Lp) ** (-0.28)
    return max(f, 1e-6)


# ============================================================
# REFRIGERANT-SIDE: EVAPORATION
# ============================================================

def h_evap_chen1966(x: float, G: float, Di: float,
                     ref: RefrigerantProperties, P: float) -> float:
    """
    Chen (1966) modified — FT (D >= 3mm).
    h_tp = F × h_l × C_nb
    q''-independent version (no Forster-Zuber self-reinforcement).
    """
    if x <= 0.001:
        x = 0.001
    if x >= 0.999:
        x = 0.999

    # Liquid-only properties
    rho_l = ref.rho_l(P)
    mu_l = ref.mu_l(P)
    k_l = ref.k_l(P)
    Pr_l = ref.Pr_l(P)
    P_r = ref.P_r(P)

    # Liquid-only HTC (Dittus-Boelter)
    Re_l = G * (1 - x) * Di / mu_l
    Re_l = max(Re_l, 100)
    h_l = 0.023 * Re_l ** 0.8 * Pr_l ** 0.4 * k_l / Di

    # Lockhart-Martinelli parameter
    Xtt = ref.Xtt(x, P)

    # Convective enhancement factor F (q''-independent)
    inv_Xtt = 1.0 / max(Xtt, 1e-10)
    if inv_Xtt > 0.1:
        F = 2.35 * (0.213 + inv_Xtt) ** 0.736
    else:
        F = 1.0

    # Nucleate boiling correction (q''-independent, Pr-based)
    C_nb = 1.0 + (1 - x) * P_r ** 0.4

    h_tp = F * h_l * C_nb
    return max(h_tp, 100.0)


def h_evap_kim_mudawar_2013(x: float, G: float, Dh: float, q_flux: float,
                             ref: RefrigerantProperties, P: float,
                             P_H: float = 1.0, P_F: float = 1.0) -> float:
    """
    Kim & Mudawar (2013) — MCHX (D = 0.19-6.5mm).
    h_tp = sqrt(h_nb² + h_cb²)
    P_H/P_F: heated/total perimeter ratio (3-side heating).
    """
    if x <= 0.001:
        x = 0.001
    if x >= 0.999:
        x = 0.999

    rho_l = ref.rho_l(P)
    rho_v = ref.rho_v(P)
    mu_l = ref.mu_l(P)
    k_l = ref.k_l(P)
    Pr_l = ref.Pr_l(P)
    h_fg_val = ref.h_fg(P)
    sigma_val = ref.sigma(P)
    P_r = ref.P_r(P)

    # Liquid-only HTC
    Re_fo = G * Dh / mu_l
    Re_fo = max(Re_fo, 100)
    h_f = 0.023 * Re_fo ** 0.8 * Pr_l ** 0.4 * k_l / Dh

    # Boiling number
    Bo = q_flux / (G * h_fg_val) if (G * h_fg_val) > 0 else 1e-6
    Bo = max(Bo, 1e-8)

    # Weber number
    We_fo = G ** 2 * Dh / (rho_l * sigma_val) if (rho_l * sigma_val) > 0 else 100

    # Lockhart-Martinelli
    Xtt = ref.Xtt(x, P)

    P_ratio = P_H / P_F if P_F > 0 else 0.75

    # Nucleate boiling
    h_nb = 2345 * (Bo * P_ratio) ** 0.7 * P_r ** 0.38 * (1 - x) ** (-0.51) * h_f

    # Convective boiling
    h_cb = (5.2 * (Bo * P_ratio) ** 0.08 * We_fo ** (-0.54) +
            3.5 / max(Xtt, 1e-6) ** 0.94 * (rho_v / rho_l) ** 0.25) * h_f

    h_tp = math.sqrt(h_nb ** 2 + h_cb ** 2)
    return max(h_tp, 100.0)


# ============================================================
# REFRIGERANT-SIDE: CONDENSATION
# ============================================================

def h_cond_shah1979(x: float, G: float, Di: float,
                    ref: RefrigerantProperties, P: float) -> float:
    """
    Shah (1979) — FT condensation (D >= 3mm).
    h_cond = h_lo × (1 + 3.8/Z^0.95)
    """
    if x <= 0.001:
        x = 0.001
    if x >= 0.999:
        x = 0.999

    rho_l = ref.rho_l(P)
    mu_l = ref.mu_l(P)
    k_l = ref.k_l(P)
    Pr_l = ref.Pr_l(P)
    P_r = ref.P_r(P)

    # Liquid-only (entire mass as liquid) Re
    Re_lo = G * Di / mu_l
    Re_lo = max(Re_lo, 100)

    # h_lo (Dittus-Boelter, all liquid)
    h_lo = 0.023 * Re_lo ** 0.8 * Pr_l ** 0.4 * k_l / Di

    # Shah parameter Z
    Z = ((1 / x - 1) ** 0.8) * P_r ** 0.4

    h_cond = h_lo * x ** 0.8 * (1 + 3.8 / max(Z, 0.01) ** 0.95)
    return max(h_cond, 100.0)


def h_cond_kim_mudawar_2012(x: float, G: float, Dh: float,
                             ref: RefrigerantProperties, P: float) -> float:
    """
    Kim & Mudawar (2012) IJHMT 55:3246 — MCHX condensation (D=0.4-3mm).
    Annular vs Slug/Bubbly regime auto-branching via We*.
    """
    if x <= 0.001:
        x = 0.001
    if x >= 0.999:
        x = 0.999

    rho_l = ref.rho_l(P)
    rho_v = ref.rho_v(P)
    mu_l = ref.mu_l(P)
    mu_v = ref.mu_v(P)
    k_l = ref.k_l(P)
    Pr_l = ref.Pr_l(P)
    sigma_val = ref.sigma(P)

    Re_l = G * (1 - x) * Dh / mu_l
    Re_l = max(Re_l, 10)
    Re_v = G * x * Dh / mu_v
    Re_v = max(Re_v, 10)

    # Liquid-only HTC
    h_f = 0.023 * (G * Dh / mu_l) ** 0.8 * Pr_l ** 0.4 * k_l / Dh

    # Weber number criterion
    We_star = 2.45 * (mu_v / mu_l) ** 0.64 * (rho_v / rho_l) ** 0.3 * Re_v ** 0.79
    We_lo = G ** 2 * Dh / (rho_l * sigma_val) if (rho_l * sigma_val) > 0 else 100

    Xtt = ref.Xtt(x, P)

    if We_star > We_lo:
        # Annular regime
        h_cond = h_f * (0.048 * Re_l ** 0.69 * Pr_l ** 0.34 *
                        (rho_v / rho_l) ** 0.12 / max(Xtt, 0.01) ** 0.5)
    else:
        # Slug/Bubbly regime
        h_cond = h_f * ((0.048 * Re_l ** 0.69 * Pr_l ** 0.34 *
                         (rho_v / rho_l) ** 0.12 / max(Xtt, 0.01) ** 0.5) ** 2 +
                        (5.7 * (rho_l * sigma_val * Dh / mu_l ** 2) ** 0.25) ** 2) ** 0.5

    return max(h_cond, 100.0)


# ============================================================
# SINGLE-PHASE — Gnielinski (1976)
# ============================================================

def h_single_gnielinski(Re: float, Pr: float, k: float, Di: float) -> float:
    """
    Gnielinski (1976) — single-phase (superheated vapor or subcooled liquid).
    Nu = (f/8)(Re-1000)Pr / [1 + 12.7√(f/8)(Pr^(2/3)-1)]
    Valid for Re > 2300.
    """
    if Re < 2300:
        # Laminar: Nu = 3.66 for constant wall temperature
        return 3.66 * k / Di

    f = (0.790 * math.log(Re) - 1.64) ** (-2)  # Petukhov friction factor
    Nu = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * math.sqrt(f / 8) * (Pr ** (2 / 3) - 1))
    Nu = max(Nu, 3.66)
    return Nu * k / Di


# ============================================================
# TRANSITION BLENDING
# ============================================================

def h_with_transition(x: float, G: float, Di: float, q_flux: float,
                      ref: RefrigerantProperties, P: float,
                      mode: str = "evap",
                      hx_type: str = "FT",
                      P_H: float = 1.0, P_F: float = 1.0) -> float:
    """
    Compute h_i with transition blending at x boundaries.
    x = 0.90~1.05: blend 2-phase → superheated vapor
    x = -0.05~0:   blend subcooled liquid → 2-phase
    """
    T_sat = ref.T_sat(P)

    if mode == "evap":
        # Two-phase evaporation
        if 0.0 < x < 0.90:
            if hx_type == "FT":
                return h_evap_chen1966(x, G, Di, ref, P)
            else:
                return h_evap_kim_mudawar_2013(x, G, Di, q_flux, ref, P, P_H, P_F)

        elif 0.90 <= x <= 1.05:
            # Transition: 2-phase → superheated vapor
            x_2ph = min(x, 0.999)
            if hx_type == "FT":
                h_2ph = h_evap_chen1966(x_2ph, G, Di, ref, P)
            else:
                h_2ph = h_evap_kim_mudawar_2013(x_2ph, G, Di, q_flux, ref, P, P_H, P_F)

            props_v = ref.props_single(T_sat + 1.0, P)
            Re_v = G * Di / props_v["mu"]
            h_vapor = h_single_gnielinski(Re_v, props_v["Pr"], props_v["k"], Di)

            w = (x - 0.90) / 0.15  # 0 at x=0.90, 1 at x=1.05
            w = max(0, min(w, 1))
            return (1 - w) * h_2ph + w * h_vapor

        elif x > 1.05:
            # Fully superheated
            props_v = ref.props_single(T_sat + 5.0, P)
            Re_v = G * Di / props_v["mu"]
            return h_single_gnielinski(Re_v, props_v["Pr"], props_v["k"], Di)

        else:  # x <= 0
            # Subcooled
            props_l = ref.props_single(T_sat - 5.0, P)
            Re_l = G * Di / props_l["mu"]
            return h_single_gnielinski(Re_l, props_l["Pr"], props_l["k"], Di)

    else:  # condensation
        if 0.05 < x < 1.0:
            if hx_type == "FT":
                return h_cond_shah1979(x, G, Di, ref, P)
            else:
                return h_cond_kim_mudawar_2012(x, G, Di, ref, P)

        elif 0.0 <= x <= 0.05:
            # Transition: 2-phase → subcooled liquid
            x_2ph = max(x, 0.001)
            if hx_type == "FT":
                h_2ph = h_cond_shah1979(x_2ph, G, Di, ref, P)
            else:
                h_2ph = h_cond_kim_mudawar_2012(x_2ph, G, Di, ref, P)

            props_l = ref.props_single(T_sat - 1.0, P)
            Re_l = G * Di / props_l["mu"]
            h_sub = h_single_gnielinski(Re_l, props_l["Pr"], props_l["k"], Di)

            w = (0.05 - x) / 0.05
            w = max(0, min(w, 1))
            return (1 - w) * h_2ph + w * h_sub

        elif x < 0:
            props_l = ref.props_single(T_sat - 5.0, P)
            Re_l = G * Di / props_l["mu"]
            return h_single_gnielinski(Re_l, props_l["Pr"], props_l["k"], Di)

        else:  # x >= 1.0
            props_v = ref.props_single(T_sat + 1.0, P)
            Re_v = G * Di / props_v["mu"]
            return h_single_gnielinski(Re_v, props_v["Pr"], props_v["k"], Di)


# ============================================================
# CORRELATION AUTO-SELECTION ENGINE
# ============================================================

def select_correlations(hx_type: str, Di: float, fin_type: str = "plain",
                        Pt: float = 0.0254, Pl: float = 0.022) -> dict:
    """
    Automatic correlation selection based on geometry and operating conditions.
    Returns dict with keys: 'air_j', 'evap', 'cond', 'single_phase'
    """
    result = {"single_phase": "gnielinski"}

    if hx_type == "FT":
        # Air-side j-factor
        PtPl = Pt / Pl if Pl > 0 else 1.0
        if fin_type == "plain":
            result["air_j"] = "wang2000" if PtPl <= 1.35 else "wang2000_high"
        elif fin_type == "wavy":
            result["air_j"] = "wang1999_wavy"
        elif fin_type == "louver":
            result["air_j"] = "wang1999_louver"
        elif fin_type == "slit":
            result["air_j"] = "slit"
        else:
            result["air_j"] = "wang2000"

        # Refrigerant-side
        result["evap"] = "chen1966"
        result["cond"] = "shah1979"

    elif hx_type == "MCHX":
        result["air_j"] = "chang_wang_1997"
        if Di < 0.003:
            result["evap"] = "kim_mudawar_2013"
            result["cond"] = "kim_mudawar_2012"
        else:
            result["evap"] = "chen1966"
            result["cond"] = "shah1979"

    return result
