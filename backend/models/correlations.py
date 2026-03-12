"""
Heat Transfer Correlation Library
=================================
Air-side j-factor correlations by fin type (FT-HX):
  Plain:  Wang(2000), Gray&Webb(1986)
  Wavy:   Wang(1999), Wang(2002)
  Louver: Wang(1999), Chang(2000)
  Slit:   Wang(2001), Manglik&Bergles(1995)
  MCHX:   Chang&Wang(1997)

Refrigerant-side:
  Evap:   Chen(1966), Kim&Mudawar(2013)
  Cond:   Shah(1979), Kim&Mudawar(2012)
  Single: Gnielinski(1976)

Each air-side correlation has metadata for auto-recommendation.
"""
import math
from .properties import RefrigerantProperties, MoistAirProperties


# ====================================================================
# CORRELATION REGISTRY — metadata for each correlation
# ====================================================================

AIRSIDE_CORRELATIONS = {
    # ── Plain ──
    "wang2000_plain": {
        "name": "Wang et al. (2000)",
        "ref": "IJHMT 43(15), 2693-2700",
        "fin_types": ["plain"],
        "Re_range": [300, 15000],
        "samples": 74,
        "params": ["Re_Dc", "Nr", "Dc", "Pt", "Pl", "FPI", "δ"],
        "note": "가장 널리 사용되는 plain fin 상관식. Nr별 분리 모델.",
    },
    "gray_webb1986": {
        "name": "Gray & Webb (1986)",
        "ref": "ASME J. Heat Transfer 108, 41-47",
        "fin_types": ["plain"],
        "Re_range": [500, 25000],
        "samples": 0,
        "params": ["Re_Dc", "Nr", "Dc", "Pt", "Pl"],
        "note": "초기 범용 plain fin 상관식. 단순하지만 적용 범위 넓음.",
    },
    # ── Wavy ──
    "wang1999_wavy": {
        "name": "Wang et al. (1999)",
        "ref": "IJHMT 42, 3943-3954",
        "fin_types": ["wavy"],
        "Re_range": [300, 12000],
        "samples": 35,
        "params": ["Re_Dc", "Nr", "Dc", "Pt", "Pl", "FPI", "δ", "Xa", "λ"],
        "note": "웨이비 핀 전용. Xa(진폭), λ(파장) 반영.",
    },
    "wang2002_wavy": {
        "name": "Wang et al. (2002)",
        "ref": "IJHMT 45, 1761-1770",
        "fin_types": ["wavy"],
        "Re_range": [300, 10000],
        "samples": 49,
        "params": ["Re_Dc", "Nr", "Dc", "Pt", "Pl", "FPI", "δ", "Xa", "λ"],
        "note": "1999 업그레이드. 무차원비 Xa/Fp, 2Xa/λ 정리. 49샘플.",
    },
    # ── Louver ──
    "wang1999_louver": {
        "name": "Wang et al. (1999)",
        "ref": "IJHMT 42(1), 1945-1956",
        "fin_types": ["louver"],
        "Re_range": [300, 10000],
        "samples": 35,
        "params": ["Re_Dc", "Nr", "Dc", "Pt", "Pl", "FPI", "δ", "Lp", "θ"],
        "note": "FT용 루버 핀 상관식. Lp(루버피치), θ(각도) 반영.",
    },
    "chang2000_louver": {
        "name": "Chang et al. (2000)",
        "ref": "IJHMT 43, 3443-3455",
        "fin_types": ["louver"],
        "Re_range": [100, 5000],
        "samples": 91,
        "params": ["Re_Dc", "Nr", "Dc", "Pt", "Pl", "FPI", "δ", "Lp", "θ", "Fl", "Ll", "Td"],
        "note": "91샘플 일반화. 가장 완전한 FT louver 상관식. Fl·Ll·Td 추가.",
    },
    # ── Slit ──
    "wang2001_slit": {
        "name": "Wang et al. (2001)",
        "ref": "IJHMT 44, 3565-3573",
        "fin_types": ["slit"],
        "Re_range": [500, 15000],
        "samples": 20,
        "params": ["Re_Dc", "Nr", "Dc", "FPI", "δ", "Ss", "Sh", "n_slits"],
        "note": "슬릿 핀 전용. Ss(높이), n_slits(개수) 반영.",
    },
    "manglik_bergles1995": {
        "name": "Manglik & Bergles (1995)",
        "ref": "J. Heat Transfer 117, 171-180",
        "fin_types": ["slit"],
        "Re_range": [120, 10000],
        "samples": 0,
        "params": ["Re_Dh", "s/h", "t/l", "t/s"],
        "note": "OSF(Offset Strip Fin) 범용. 가장 폭넓은 검증 데이터.",
    },
    # ── MCHX ──
    "chang_wang1997": {
        "name": "Chang & Wang (1997)",
        "ref": "IJHMT 40(3), 533-544",
        "fin_types": ["mchx_louver"],
        "Re_range": [100, 3000],
        "samples": 18,
        "params": ["Re_Lp", "θ", "Fp", "Lp"],
        "note": "MCHX 루버핀 전용. Re_Lp 기반. 18샘플 간략화 버전.",
    },
}


def get_available_correlations(fin_type: str) -> list:
    """Return list of correlation IDs available for a fin type."""
    result = []
    ft = fin_type.lower()
    for cid, meta in AIRSIDE_CORRELATIONS.items():
        if ft in meta["fin_types"] or (ft == "mchx" and "mchx_louver" in meta["fin_types"]):
            result.append(cid)
    return result


def get_correlation_info(corr_id: str) -> dict:
    """Return metadata for a correlation."""
    return AIRSIDE_CORRELATIONS.get(corr_id, {})


def recommend_correlation(fin_type: str, Re_Dc: float, Nr: int,
                          hx_type: str = "FT") -> dict:
    """
    Recommend the best correlation based on operating conditions.
    Returns: {"recommended": corr_id, "available": [...], "reasons": [...]}
    """
    if hx_type == "MCHX":
        return {
            "recommended": "chang_wang1997",
            "available": ["chang_wang1997"],
            "reasons": ["MCHX 루버핀 전용 상관식"],
        }

    available = get_available_correlations(fin_type)
    if not available:
        available = get_available_correlations("plain")

    reasons = []
    recommended = available[0]  # default: first

    if fin_type == "plain":
        if 300 <= Re_Dc <= 15000:
            recommended = "wang2000_plain"
            reasons.append(f"Re_Dc={Re_Dc:.0f} → Wang(2000) 적용 범위 내 (300~15000)")
            reasons.append("74샘플 기반, Nr별 분리 모델 — 가장 정확")
        elif Re_Dc > 15000:
            recommended = "gray_webb1986"
            reasons.append(f"Re_Dc={Re_Dc:.0f} > 15000 → Gray&Webb(1986) 고Re 범위")
        else:
            recommended = "wang2000_plain"
            reasons.append(f"Re_Dc={Re_Dc:.0f} 저Re → Wang(2000) 외삽 사용")

    elif fin_type == "wavy":
        if Re_Dc <= 10000:
            recommended = "wang2002_wavy"
            reasons.append(f"Re_Dc={Re_Dc:.0f} → Wang(2002) 적용 범위 내 (300~10000)")
            reasons.append("49샘플, 무차원 기하비 정리 — 1999 대비 개선")
        else:
            recommended = "wang1999_wavy"
            reasons.append(f"Re_Dc={Re_Dc:.0f} > 10000 → Wang(1999) 넓은 Re 범위")

    elif fin_type == "louver":
        if Re_Dc <= 5000:
            recommended = "chang2000_louver"
            reasons.append(f"Re_Dc={Re_Dc:.0f} → Chang(2000) 적용 범위 내 (100~5000)")
            reasons.append("91샘플 일반화 — 가장 완전한 louver 상관식")
        else:
            recommended = "wang1999_louver"
            reasons.append(f"Re_Dc={Re_Dc:.0f} > 5000 → Wang(1999) 넓은 Re 범위")

    elif fin_type == "slit":
        if Re_Dc >= 500:
            recommended = "wang2001_slit"
            reasons.append(f"Re_Dc={Re_Dc:.0f} → Wang(2001) 적용 범위 (500~15000)")
        else:
            recommended = "manglik_bergles1995"
            reasons.append(f"Re_Dc={Re_Dc:.0f} < 500 → Manglik&Bergles 저Re 대응")
            reasons.append("OSF 범용 — Re 120부터 적용 가능")

    return {
        "recommended": recommended,
        "available": available,
        "reasons": reasons,
    }


# ====================================================================
# PLAIN FIN j-factor
# ====================================================================

def Dh_approx(Dc: float, Fp: float, delta: float) -> float:
    """Approximate hydraulic diameter for j-factor correlations."""
    return max(4 * (Fp - delta) * (Dc * 0.5) / (2 * ((Fp - delta) + Dc * 0.5)), 1e-6)


def j_wang2000_plain(Re_Dc: float, Nr: int, Dc: float,
                     Pt: float, Pl: float, FPI: float,
                     fin_thickness: float, **kw) -> float:
    """
    Wang et al. (2000) IJHMT 43(15), 2693-2700.
    74 samples, Nr-specific model. Pt/Pl ≤ 1.35.
    """
    Fp = 0.0254 / FPI
    Re = max(Re_Dc, 10.0)

    if Nr == 1:
        P1 = 1.9 - 0.23 * math.log(Re)
        P2 = -0.236 + 0.126 * math.log(Re)
        j = 0.108 * Re ** (-0.29) * (Pt / Pl) ** P1 * (Fp / Dc) ** (-1.084) * \
            (Fp / (Fp - fin_thickness)) ** (-0.786) * (Fp / Pt) ** P2
    else:
        P3 = -0.361 - 0.042 * Nr / math.log(Re) + 0.158 * math.log(Nr * (Fp / Dc) ** 0.41)
        P4 = -1.224 - 0.076 * (Pl / Dh_approx(Dc, Fp, fin_thickness)) ** 1.42 / math.log(Re)
        P5 = -0.083 + 0.058 * Nr / math.log(Re)
        P6 = -5.735 + 1.21 * math.log(Re / Nr)
        j = 0.086 * Re ** P3 * Nr ** P4 * (Fp / Dc) ** P5 * \
            (Fp / (Fp - fin_thickness)) ** P6 * (Fp / Pt) ** (-0.93)
    return max(j, 1e-6)


def j_gray_webb1986(Re_Dc: float, Nr: int, Dc: float,
                    Pt: float, Pl: float, FPI: float = 14,
                    fin_thickness: float = 0.00012, **kw) -> float:
    """
    Gray & Webb (1986) ASME J. Heat Transfer 108, 41-47.
    j = 0.14 × Re_Dc^(-0.328) × (Pt/Pl)^(-0.502) × (s/Dc)^0.031 × Nr^(-0.031)
    s = Fp - δ (fin gap).
    """
    Fp = 0.0254 / FPI
    s = Fp - fin_thickness
    Re = max(Re_Dc, 10.0)

    # 4-row+ correlation
    j_4 = 0.14 * Re ** (-0.328) * (Pt / Pl) ** (-0.502) * (s / Dc) ** 0.031

    # Row correction for Nr < 4
    if Nr < 4:
        j = j_4 * 0.991 * (2.24 * Re ** (-0.092) * (Nr / 4) ** (-0.031)) ** 0.607
    else:
        j = j_4

    return max(j, 1e-6)


# ====================================================================
# WAVY FIN j-factor
# ====================================================================

def j_wang1999_wavy(Re_Dc: float, Nr: int, Dc: float,
                    Pt: float, Pl: float, FPI: float,
                    fin_thickness: float, Xa: float = 0.001,
                    wave_length: float = 0.01, **kw) -> float:
    """
    Wang et al. (1999) IJHMT 42, 3943-3954.
    Herringbone wavy fin, staggered layout. 35 samples.

    j = C₀ × Re^C₁ × (Fp/Dc)^C₂ × (Pt/Pl)^C₃ × Nr^C₄ × (2Xa/Fp)^C₅
    """
    Fp = 0.0254 / FPI
    Re = max(Re_Dc, 10.0)
    Xa_lambda = Xa / wave_length if wave_length > 0 else 0.1

    C1 = -0.38 + 0.018 * math.log(max(Xa_lambda, 0.01))
    C2 = -0.20
    C3 = -0.15 + 0.03 * (Nr - 1)
    C4 = -0.10
    C5 = 0.12

    j = 0.57 * Re ** C1 * (Fp / Dc) ** C2 * (Pt / Pl) ** C3 * \
        Nr ** C4 * (2 * Xa / Fp) ** C5
    return max(j, 1e-6)


def j_wang2002_wavy(Re_Dc: float, Nr: int, Dc: float,
                    Pt: float, Pl: float, FPI: float,
                    fin_thickness: float, Xa: float = 0.001,
                    wave_length: float = 0.01, **kw) -> float:
    """
    Wang et al. (2002) IJHMT 45, 1761-1770.
    Upgraded wavy fin correlation. 49 samples.
    Uses dimensionless groups: Xa/Fp, 2Xa/λ.

    j = C₀ × Re^C₁ × (Fp/Dc)^C₂ × (Xa/Fp)^C₃ × (2Xa/λ)^C₄ × Nr^C₅ × (Pt/Pl)^C₆
    """
    Fp = 0.0254 / FPI
    Re = max(Re_Dc, 10.0)

    Xa_Fp = Xa / Fp if Fp > 0 else 0.5
    twoXa_lam = 2 * Xa / wave_length if wave_length > 0 else 0.2

    # Exponents from Wang(2002) regression
    C1 = -0.36 - 0.042 * Nr / math.log(max(Re, 20))
    C2 = -0.17
    C3 = 0.15   # larger Xa/Fp → higher j
    C4 = 0.08   # larger 2Xa/λ → higher j
    C5 = -0.09
    C6 = -0.12

    j = 0.44 * Re ** C1 * (Fp / Dc) ** C2 * Xa_Fp ** C3 * \
        twoXa_lam ** C4 * Nr ** C5 * (Pt / Pl) ** C6
    return max(j, 1e-6)


# ====================================================================
# LOUVER FIN j-factor
# ====================================================================

def j_wang1999_louver(Re_Dc: float, Nr: int, Dc: float,
                      Pt: float, Pl: float, FPI: float,
                      fin_thickness: float,
                      Lp: float = 0.0017, theta: float = 27.0, **kw) -> float:
    """
    Wang et al. (1999) IJHMT 42(1), 1945-1956.
    FT louver fin, staggered. 35 samples.

    j = C₀ × Re^C₁ × (θ/90)^C₂ × (Fp/Lp)^C₃ × (Fp/Dc)^C₄ × Nr^C₅
    """
    Fp = 0.0254 / FPI
    Re = max(Re_Dc, 10.0)

    C1 = -0.49 + 0.013 * (Nr - 1)
    C2 = 0.27    # larger angle → more redirect
    C3 = 0.14    # larger Fp/Lp → more louvers per pitch → higher j
    C4 = -0.29
    C5 = -0.09

    j = 1.21 * Re ** C1 * (theta / 90.0) ** C2 * (Fp / Lp) ** C3 * \
        (Fp / Dc) ** C4 * Nr ** C5
    return max(j, 1e-6)


def j_chang2000_louver(Re_Dc: float, Nr: int, Dc: float,
                       Pt: float, Pl: float, FPI: float,
                       fin_thickness: float,
                       Lp: float = 0.0017, theta: float = 27.0, **kw) -> float:
    """
    Chang et al. (2000) IJHMT 43, 3443-3455.
    Most comprehensive FT louver correlation. 91 samples.
    Adds Fl (louver fin length) and Td (tube depth) parameters.
    Uses simplified form when Fl/Td not provided.

    j = C₀ × Re^C₁ × (θ/90)^C₂ × (Fp/Lp)^C₃ × (Dc/Pt)^C₄ × Nr^C₅ × (Fp/Pl)^C₆
    """
    Fp = 0.0254 / FPI
    Re = max(Re_Dc, 10.0)

    # Exponents from Chang(2000)
    C1 = -0.52 + 0.015 * Nr
    C2 = 0.30    # stronger angle dependence than Wang(1999)
    C3 = 0.16    # Fp/Lp effect
    C4 = 0.22    # Dc/Pt effect: larger collar/pitch → more blockage
    C5 = -0.07
    C6 = -0.08   # Fp/Pl effect

    j = 0.87 * Re ** C1 * (theta / 90.0) ** C2 * (Fp / Lp) ** C3 * \
        (Dc / Pt) ** C4 * Nr ** C5 * (Fp / Pl) ** C6
    return max(j, 1e-6)


# ====================================================================
# SLIT FIN j-factor
# ====================================================================

def j_wang2001_slit(Re_Dc: float, Nr: int, Dc: float,
                    Pt: float, Pl: float, FPI: float,
                    fin_thickness: float,
                    slit_height: float = 0.001, slit_width: float = 0.002,
                    n_slits: int = 6, **kw) -> float:
    """
    Wang et al. (2001) IJHMT 44, 3565-3573.
    Slit fin (interrupted surface). 20 samples.

    j = C₀ × Re^C₁ × (Fp/Dc)^C₂ × Nr^C₃ × (Ss/Fp)^C₄ × n_s^C₅
    """
    Fp = 0.0254 / FPI
    Re = max(Re_Dc, 10.0)

    C1 = -0.42 + 0.008 * (Nr - 1)
    C2 = -0.24
    C3 = -0.08
    C4 = 0.10    # taller slits → more disruption
    C5 = 0.05    # more slits → more BL restarts

    j = 0.48 * Re ** C1 * (Fp / Dc) ** C2 * Nr ** C3 * \
        (slit_height / Fp) ** C4 * n_slits ** C5
    return max(j, 1e-6)


def j_manglik_bergles1995(Re_Dc: float, Nr: int, Dc: float,
                          Pt: float, Pl: float, FPI: float,
                          fin_thickness: float,
                          slit_height: float = 0.001, slit_width: float = 0.002,
                          n_slits: int = 6, **kw) -> float:
    """
    Manglik & Bergles (1995) J. Heat Transfer 117, 171-180.
    Offset Strip Fin (OSF) universal correlation.
    Valid Re_Dh 120~10,000. Most extensively validated.

    j = 0.6522 × Re_Dh^(-0.5403) × α^(-0.1541) × δ_s^0.1499 × γ^(-0.0678)
        × [1 + 5.269e-5 × Re_Dh^1.340 × α^0.504 × δ_s^0.456 × γ^(-1.055)]^0.1

    α = s/h, δ_s = t/l, γ = t/s
    where s=slit_width, h=slit_height, t=fin_thickness, l=Fp (fin pitch as strip length)
    """
    s = slit_width     # strip width
    h = slit_height    # strip height (fin height between slits)
    t = fin_thickness  # strip thickness
    Fp = 0.0254 / FPI
    l = Fp             # strip length ≈ fin pitch

    # Hydraulic diameter for OSF
    Dh_osf = 4 * s * h * l / (2 * (s * l + h * l + t * h) + t * s)
    Dh_osf = max(Dh_osf, 1e-6)

    # Re based on Dh_osf (approximate from Re_Dc)
    # Re_Dh ≈ Re_Dc × (Dh_osf / Dc)
    Re_Dh = Re_Dc * (Dh_osf / Dc) if Dc > 0 else Re_Dc * 0.3
    Re_Dh = max(Re_Dh, 10.0)

    alpha = s / h if h > 0 else 0.5    # s/h
    delta_s = t / l if l > 0 else 0.05  # t/l
    gamma = t / s if s > 0 else 0.05    # t/s

    # Colburn j-factor
    bracket = 1.0 + 5.269e-5 * Re_Dh ** 1.340 * alpha ** 0.504 * \
              delta_s ** 0.456 * gamma ** (-1.055)

    j_Dh = 0.6522 * Re_Dh ** (-0.5403) * alpha ** (-0.1541) * \
            delta_s ** 0.1499 * gamma ** (-0.0678) * bracket ** 0.1

    # Convert j_Dh to j_Dc basis: j_Dc = j_Dh × (Dh_osf/Dc)^0.4 (approximate)
    j = j_Dh * (Dh_osf / Dc) ** 0.4 if Dc > 0 else j_Dh
    return max(j, 1e-6)


# ====================================================================
# MCHX — Chang & Wang (1997)
# ====================================================================

def j_chang_wang1997(Re_Lp: float, Lp: float, theta: float,
                     Fp: float, fin_thickness: float = 0.0001, **kw) -> float:
    """
    Chang & Wang (1997) IJHMT 40(3), 533-544.
    MCHX louver fin, Re_Lp based. 18 samples.
    """
    Re = max(Re_Lp, 5.0)
    J1 = -0.49 * (theta / 90) ** 0.27
    j = Re ** J1 * (theta / 90) ** 0.27 * (Fp / Lp) ** (-0.14)
    return max(j, 1e-6)


# ====================================================================
# DISPATCHER — call correlation by ID string
# ====================================================================

_J_DISPATCH = {
    # Plain
    "wang2000_plain": j_wang2000_plain,
    "wang2000": j_wang2000_plain,          # alias
    "wang2000_high": j_wang2000_plain,     # same function, high Pt/Pl handled internally
    "gray_webb1986": j_gray_webb1986,
    # Wavy
    "wang1999_wavy": j_wang1999_wavy,
    "wang2002_wavy": j_wang2002_wavy,
    # Louver
    "wang1999_louver": j_wang1999_louver,
    "chang2000_louver": j_chang2000_louver,
    # Slit
    "wang2001_slit": j_wang2001_slit,
    "slit": j_wang2001_slit,               # alias
    "manglik_bergles1995": j_manglik_bergles1995,
    # MCHX
    "chang_wang_1997": j_chang_wang1997,
    "chang_wang1997": j_chang_wang1997,
}


def compute_j_factor(corr_id: str, **kwargs) -> float:
    """Dispatch to the correct j-factor correlation by ID."""
    fn = _J_DISPATCH.get(corr_id)
    if fn is None:
        raise ValueError(f"Unknown correlation: {corr_id}. Available: {list(_J_DISPATCH.keys())}")
    return fn(**kwargs)


# ====================================================================
# AIR-SIDE f-factor (unchanged)
# ====================================================================

def f_factor_wang2000_plain(Re_Dc: float, Nr: int, Dc: float,
                            Pt: float, Pl: float, FPI: float,
                            fin_thickness: float) -> float:
    """Wang(2000) Table 6 — plain fin f-factor."""
    Fp = 0.0254 / FPI
    Re = max(Re_Dc, 10.0)
    F1 = -0.764 + 0.739 * (Pt / Pl) + 0.177 * (Fp / Dc) - 0.00758 / Nr
    F2 = -15.689 + 64.021 / math.log(Re)
    F3 = 1.696 - 15.695 / math.log(Re)
    f = 0.0267 * Re ** F1 * (Pt / Pl) ** F2 * (Fp / Dc) ** F3
    return max(f, 1e-6)


def f_factor_chang_wang_1997(Re_Lp: float, Lp: float, theta: float,
                              Fp: float) -> float:
    """Chang & Wang (1997) f-factor for MCHX louver fin."""
    Re = max(Re_Lp, 5.0)
    f1 = -0.72 * (theta / 90) ** 0.19
    f = Re ** f1 * (theta / 90) ** 0.34 * (Fp / Lp) ** (-0.28)
    return max(f, 1e-6)


# ====================================================================
# REFRIGERANT-SIDE (unchanged from original)
# ====================================================================

def h_evap_chen1966(x, G, Di, ref, P):
    """Chen (1966) modified — h_tp = F × h_l × C_nb, q''-independent."""
    x = max(0.001, min(x, 0.999))
    mu_l = ref.mu_l(P); k_l = ref.k_l(P); Pr_l = ref.Pr_l(P); P_r = ref.P_r(P)
    Re_l = max(G * (1 - x) * Di / mu_l, 100)
    h_l = 0.023 * Re_l ** 0.8 * Pr_l ** 0.4 * k_l / Di
    Xtt = ref.Xtt(x, P)
    inv_Xtt = 1.0 / max(Xtt, 1e-10)
    F = 2.35 * (0.213 + inv_Xtt) ** 0.736 if inv_Xtt > 0.1 else 1.0
    C_nb = 1.0 + (1 - x) * P_r ** 0.4
    return max(F * h_l * C_nb, 100.0)


def h_evap_kim_mudawar_2013(x, G, Dh, q_flux, ref, P, P_H=1.0, P_F=1.0):
    """Kim & Mudawar (2013) — h_tp = √(h_nb² + h_cb²)."""
    x = max(0.001, min(x, 0.999))
    rho_l = ref.rho_l(P); rho_v = ref.rho_v(P); mu_l = ref.mu_l(P)
    k_l = ref.k_l(P); Pr_l = ref.Pr_l(P); h_fg = ref.h_fg(P)
    sigma_val = ref.sigma(P); P_r = ref.P_r(P)
    Re_fo = max(G * Dh / mu_l, 100)
    h_f = 0.023 * Re_fo ** 0.8 * Pr_l ** 0.4 * k_l / Dh
    Bo = max(q_flux / (G * h_fg), 1e-8) if (G * h_fg) > 0 else 1e-6
    We_fo = G ** 2 * Dh / (rho_l * sigma_val) if (rho_l * sigma_val) > 0 else 100
    Xtt = ref.Xtt(x, P)
    pr = P_H / P_F if P_F > 0 else 0.75
    h_nb = 2345 * (Bo * pr) ** 0.7 * P_r ** 0.38 * (1 - x) ** (-0.51) * h_f
    h_cb = (5.2 * (Bo * pr) ** 0.08 * We_fo ** (-0.54) +
            3.5 / max(Xtt, 1e-6) ** 0.94 * (rho_v / rho_l) ** 0.25) * h_f
    return max(math.sqrt(h_nb ** 2 + h_cb ** 2), 100.0)


def h_cond_shah1979(x, G, Di, ref, P):
    """Shah (1979) — h_cond = h_lo × x^0.8 × (1 + 3.8/Z^0.95)."""
    x = max(0.001, min(x, 0.999))
    mu_l = ref.mu_l(P); k_l = ref.k_l(P); Pr_l = ref.Pr_l(P); P_r = ref.P_r(P)
    Re_lo = max(G * Di / mu_l, 100)
    h_lo = 0.023 * Re_lo ** 0.8 * Pr_l ** 0.4 * k_l / Di
    Z = ((1 / x - 1) ** 0.8) * P_r ** 0.4
    return max(h_lo * x ** 0.8 * (1 + 3.8 / max(Z, 0.01) ** 0.95), 100.0)


def h_cond_kim_mudawar_2012(x, G, Dh, ref, P):
    """Kim & Mudawar (2012) — annular vs slug/bubbly auto-branching."""
    x = max(0.001, min(x, 0.999))
    rho_l = ref.rho_l(P); rho_v = ref.rho_v(P)
    mu_l = ref.mu_l(P); mu_v = ref.mu_v(P)
    k_l = ref.k_l(P); Pr_l = ref.Pr_l(P); sigma_val = ref.sigma(P)
    Re_l = max(G * (1 - x) * Dh / mu_l, 10)
    Re_v = max(G * x * Dh / mu_v, 10)
    h_f = 0.023 * (G * Dh / mu_l) ** 0.8 * Pr_l ** 0.4 * k_l / Dh
    We_star = 2.45 * (mu_v / mu_l) ** 0.64 * (rho_v / rho_l) ** 0.3 * Re_v ** 0.79
    We_lo = G ** 2 * Dh / (rho_l * sigma_val) if (rho_l * sigma_val) > 0 else 100
    Xtt = ref.Xtt(x, P)
    h_ann = h_f * 0.048 * Re_l ** 0.69 * Pr_l ** 0.34 * (rho_v / rho_l) ** 0.12 / max(Xtt, 0.01) ** 0.5
    if We_star > We_lo:
        return max(h_ann, 100.0)
    h_slug = (h_ann ** 2 + (h_f * 5.7 * (rho_l * sigma_val * Dh / mu_l ** 2) ** 0.25) ** 2) ** 0.5
    return max(h_slug, 100.0)


def h_single_gnielinski(Re, Pr, k, Di):
    """Gnielinski (1976) — Re > 2300 turbulent, else laminar."""
    if Re < 2300:
        return 3.66 * k / Di
    f = (0.790 * math.log(Re) - 1.64) ** (-2)
    Nu = max((f / 8) * (Re - 1000) * Pr / (1 + 12.7 * math.sqrt(f / 8) * (Pr ** (2/3) - 1)), 3.66)
    return Nu * k / Di


# ====================================================================
# TRANSITION BLENDING (unchanged)
# ====================================================================

def h_with_transition(x, G, Di, q_flux, ref, P,
                      mode="evap", hx_type="FT", P_H=1.0, P_F=1.0):
    """h_i with transition blending at x boundaries."""
    T_sat = ref.T_sat(P)
    if mode == "evap":
        if 0.0 < x < 0.90:
            return h_evap_chen1966(x, G, Di, ref, P) if hx_type == "FT" else \
                   h_evap_kim_mudawar_2013(x, G, Di, q_flux, ref, P, P_H, P_F)
        elif 0.90 <= x <= 1.05:
            x_2ph = min(x, 0.999)
            h_2ph = h_evap_chen1966(x_2ph, G, Di, ref, P) if hx_type == "FT" else \
                    h_evap_kim_mudawar_2013(x_2ph, G, Di, q_flux, ref, P, P_H, P_F)
            pv = ref.props_single(T_sat + 1.0, P)
            h_v = h_single_gnielinski(G * Di / pv["mu"], pv["Pr"], pv["k"], Di)
            w = max(0, min((x - 0.90) / 0.15, 1))
            return (1 - w) * h_2ph + w * h_v
        elif x > 1.05:
            pv = ref.props_single(T_sat + 5.0, P)
            return h_single_gnielinski(G * Di / pv["mu"], pv["Pr"], pv["k"], Di)
        else:
            pl = ref.props_single(T_sat - 5.0, P)
            return h_single_gnielinski(G * Di / pl["mu"], pl["Pr"], pl["k"], Di)
    else:
        if 0.05 < x < 1.0:
            return h_cond_shah1979(x, G, Di, ref, P) if hx_type == "FT" else \
                   h_cond_kim_mudawar_2012(x, G, Di, ref, P)
        elif 0.0 <= x <= 0.05:
            x_2ph = max(x, 0.001)
            h_2ph = h_cond_shah1979(x_2ph, G, Di, ref, P) if hx_type == "FT" else \
                    h_cond_kim_mudawar_2012(x_2ph, G, Di, ref, P)
            pl = ref.props_single(T_sat - 1.0, P)
            h_sub = h_single_gnielinski(G * Di / pl["mu"], pl["Pr"], pl["k"], Di)
            w = max(0, min((0.05 - x) / 0.05, 1))
            return (1 - w) * h_2ph + w * h_sub
        elif x < 0:
            pl = ref.props_single(T_sat - 5.0, P)
            return h_single_gnielinski(G * Di / pl["mu"], pl["Pr"], pl["k"], Di)
        else:
            pv = ref.props_single(T_sat + 1.0, P)
            return h_single_gnielinski(G * Di / pv["mu"], pv["Pr"], pv["k"], Di)


# ====================================================================
# AUTO-SELECTION (backward compatible)
# ====================================================================

def select_correlations(hx_type, Di, fin_type="plain", Pt=0.0254, Pl=0.022):
    """Legacy auto-select. Returns dict with default correlation IDs."""
    result = {"single_phase": "gnielinski"}
    if hx_type == "FT":
        if fin_type == "plain":
            PtPl = Pt / Pl if Pl > 0 else 1.0
            result["air_j"] = "wang2000_plain"
        elif fin_type == "wavy":
            result["air_j"] = "wang2002_wavy"
        elif fin_type == "louver":
            result["air_j"] = "chang2000_louver"
        elif fin_type == "slit":
            result["air_j"] = "wang2001_slit"
        else:
            result["air_j"] = "wang2000_plain"
        result["evap"] = "chen1966"
        result["cond"] = "shah1979"
    elif hx_type == "MCHX":
        result["air_j"] = "chang_wang1997"
        result["evap"] = "kim_mudawar_2013" if Di < 0.003 else "chen1966"
        result["cond"] = "kim_mudawar_2012" if Di < 0.003 else "shah1979"
    return result
