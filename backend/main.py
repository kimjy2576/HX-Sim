"""
FastAPI Backend — Heat Exchanger Simulator API
"""
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
from typing import Optional, Literal, List, Dict
from pathlib import Path
import traceback

from models.solver import SimulationInput, HXSolver
from models.geometry import FinTubeSpec, MCHXSpec
from models.properties import RefrigerantProperties
from models.correlations import (
    AIRSIDE_CORRELATIONS, get_available_correlations,
    recommend_correlation, select_correlations,
    validate_correlation, build_spec_values, build_mchx_spec_values,
    FSIDE_CORRELATIONS, get_available_f_correlations,
    recommend_f_correlation,
    REFSIDE_EVAP_CORRELATIONS, REFSIDE_COND_CORRELATIONS,
    REFSIDE_DP_CORRELATIONS,
    get_available_ref_correlations, get_available_dp_correlations,
)

app = FastAPI(
    title="HX Simulator API",
    description="열교환기 해석 모델 — FT-HX / MCHX Level 2 Tube-Segment",
    version="7.9",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ============================================================
# Pydantic Models for API
# ============================================================

class FTSpecInput(BaseModel):
    W: float = 0.5
    H: float = 0.3
    D: float = 0.08
    Do: float = 0.00952
    Di: float = 0.00822
    Pt: float = 0.0254
    Pl: float = 0.022
    Nr: int = 4
    Nt: int = 12
    layout: Literal["staggered", "inline"] = "staggered"
    FPI: float = 14.0
    fin_thickness: float = 0.00012
    fin_type: Literal["plain", "wavy", "louver", "slit"] = "plain"
    k_fin: float = 200.0
    # Wavy
    wavy_amplitude: float = 0.001
    wavy_wavelength: float = 0.01
    # Louver
    louver_pitch: float = 0.0017
    louver_angle: float = 27.0
    # Slit
    slit_height: float = 0.001
    slit_width: float = 0.002
    n_slits: int = 6
    N_seg: int = 5
    circuit_mode: str = "row_parallel"  # row_parallel, serpentine_2, serpentine_4, single, custom
    circuits: Optional[List] = None  # custom: [[row,col],[row,col],...] per circuit


class MCHXSpecInput(BaseModel):
    W: float = 0.6
    H: float = 0.4
    D: float = 0.020
    Nr: int = 1
    ch_width: float = 0.001
    ch_height: float = 0.0015
    ch_wall: float = 0.0003
    n_ports: int = 12
    tube_height: float = 0.002
    tube_pitch: float = 0.010
    louver_pitch: float = 0.0013
    louver_angle: float = 27.0
    fin_pitch: float = 0.0014
    fin_thickness: float = 0.0001
    fin_height: float = 0.008
    k_fin: float = 200.0
    N_seg: int = 5
    Nt: int = 40
    passes: Optional[List[List[int]]] = None  # baffle passes: [[tube_indices], ...]
    pass_slabs: Optional[List[int]] = None    # which slab each pass belongs to


class SimRequest(BaseModel):
    hx_type: Literal["FT", "MCHX"] = "FT"
    mode: Literal["evap", "cond"] = "evap"

    # Air
    T_air_in_C: float = Field(35.0, description="Air inlet temperature [°C]")
    RH_in: float = Field(0.50, ge=0, le=1, description="Relative humidity [-]")
    air_flow_mode: Literal["velocity", "CMM"] = "velocity"
    V_air: Optional[float] = Field(None, description="Face velocity [m/s]")
    CMM: Optional[float] = Field(None, description="Volumetric flow rate [m³/min]")

    # Refrigerant
    fluid: str = "R410A"
    T_sat_C: Optional[float] = Field(None, description="Saturation temperature [°C] (either T or P)")
    P_sat_kPa: Optional[float] = Field(None, description="Saturation pressure [kPa] (either T or P)")
    m_ref: float = Field(0.02, description="Mass flow rate [kg/s]")

    # Two-phase inlet: if True, x_in is used; if False, T_ref_in_C + P_ref_kPa define single-phase state
    two_phase_inlet: bool = Field(True, description="Is inlet two-phase?")
    x_in: float = Field(0.2, description="Inlet quality (used when two_phase_inlet=True)")
    T_ref_in_C: Optional[float] = Field(None, description="Inlet ref temp [°C] (single-phase)")

    # Flow
    flow_arrangement: Literal["counter", "parallel"] = "counter"

    # Solver settings
    max_outer: int = Field(30, description="Max outer iterations (air-ref coupling)")
    outer_tol_pct: float = Field(0.1, description="Outer convergence tolerance [%]")

    # Correlation selection (None = auto-recommend)
    air_j_corr: Optional[str] = Field(None, description="Air-side j-factor correlation ID. None = auto.")
    air_f_corr: Optional[str] = Field(None, description="Air-side f-factor correlation ID. None = auto.")
    evap_corr: Optional[str] = Field(None, description="Evaporation correlation ID. None = auto.")
    cond_corr: Optional[str] = Field(None, description="Condensation correlation ID. None = auto.")
    dp_ref_corr: Optional[str] = Field(None, description="Refrigerant dp correlation ID. None = auto.")

    # Correction factors
    cf_j: float = Field(1.0, description="Air-side j-factor correction multiplier")
    cf_f: float = Field(1.0, description="Air-side f-factor correction multiplier")
    cf_hi: float = Field(1.0, description="Refrigerant-side HTC correction multiplier")
    cf_dp_ref: float = Field(1.0, description="Refrigerant-side dp correction multiplier")

    # Geometry
    ft_spec: Optional[FTSpecInput] = None
    mchx_spec: Optional[MCHXSpecInput] = None


class SegmentOut(BaseModel):
    row: int
    tube: int
    seg: int
    Q: float
    Q_sensible: float
    Q_latent: float
    T_wall: float
    T_air_local: float
    T_ref: float
    x_ref: float
    h_i: float
    h_o: float
    eta_o: float
    dp_ref: float
    is_wet: bool
    converged: bool


class SimResponse(BaseModel):
    Q_total: float
    Q_sensible: float
    Q_latent: float
    SHR: float
    T_air_out_C: float
    W_air_out: float
    RH_out: float
    x_ref_out: float
    T_ref_out_C: float
    dp_air: float
    dp_ref: float = 0.0
    T_sat_C: float = 0.0
    P_sat_kPa: float = 0.0
    V_air: float = 0.0
    CMM: float = 0.0
    Re_Dc: float = 0.0
    row_Q: List[float]
    correlations_used: Dict
    correlation_recommendation: Dict = {}
    f_recommendation: Dict = {}
    convergence: Dict = {}
    segments: List[SegmentOut]
    error: str = ""


# ============================================================
# Endpoints
# ============================================================

STATIC_DIR = Path(__file__).parent / "static"


@app.get("/")
def root():
    """Serve frontend UI."""
    index_file = STATIC_DIR / "index.html"
    if index_file.exists():
        return FileResponse(index_file)
    return {"error": "index.html not found", "static_dir": str(STATIC_DIR),
            "files": [str(f) for f in STATIC_DIR.iterdir()] if STATIC_DIR.exists() else "dir missing"}


@app.get("/health")
def health():
    return {"status": "ok", "static_exists": (STATIC_DIR / "index.html").exists()}


@app.get("/api")
def api_info():
    return {
        "name": "HX Simulator API v7.9",
        "description": "FT-HX / MCHX Level 2 Tube-Segment Model",
        "endpoints": ["/simulate", "/refrigerants", "/docs"],
    }


@app.get("/refrigerants")
def get_refrigerants():
    return {"refrigerants": RefrigerantProperties.SUPPORTED}


@app.post("/simulate", response_model=SimResponse)
def simulate(req: SimRequest):
    """Run heat exchanger simulation."""
    try:
        # Convert °C → K
        T_air_K = req.T_air_in_C + 273.15

        # --- Resolve refrigerant state: T_sat ↔ P_sat ---
        ref_props = RefrigerantProperties(req.fluid)

        if req.T_sat_C is not None and req.P_sat_kPa is not None:
            # Both given — use T_sat as primary, P is informational
            T_sat_K = req.T_sat_C + 273.15
            P_sat = ref_props.P_sat(T_sat_K)
        elif req.T_sat_C is not None:
            T_sat_K = req.T_sat_C + 273.15
            P_sat = ref_props.P_sat(T_sat_K)
        elif req.P_sat_kPa is not None:
            P_sat = req.P_sat_kPa * 1000.0  # kPa → Pa
            T_sat_K = ref_props.T_sat(P_sat)
        else:
            # Default fallback
            T_sat_K = 280.15  # 7°C
            P_sat = ref_props.P_sat(T_sat_K)

        T_sat_C_resolved = T_sat_K - 273.15
        P_sat_kPa_resolved = P_sat / 1000.0

        # --- Determine inlet quality ---
        if req.two_phase_inlet:
            x_in = req.x_in
        else:
            # Single-phase: determine x from temperature relative to T_sat
            if req.T_ref_in_C is not None:
                T_ref_in_K = req.T_ref_in_C + 273.15
                if req.mode == "evap":
                    # Subcooled liquid entering evaporator
                    x_in = 0.0
                else:
                    # Superheated vapor entering condenser
                    x_in = 1.0
            else:
                x_in = 0.2 if req.mode == "evap" else 0.95

        # Build FT spec
        ft = None
        if req.hx_type == "FT":
            ft_in = req.ft_spec or FTSpecInput()
            ft = FinTubeSpec(
                W=ft_in.W, H=ft_in.H, D=ft_in.D,
                Do=ft_in.Do, Di=ft_in.Di,
                Pt=ft_in.Pt, Pl=ft_in.Pl,
                Nr=ft_in.Nr, Nt=ft_in.Nt,
                layout=ft_in.layout,
                FPI=ft_in.FPI, fin_thickness=ft_in.fin_thickness,
                fin_type=ft_in.fin_type, k_fin=ft_in.k_fin,
                wavy_amplitude=ft_in.wavy_amplitude, wavy_wavelength=ft_in.wavy_wavelength,
                louver_pitch=ft_in.louver_pitch, louver_angle=ft_in.louver_angle,
                slit_height=ft_in.slit_height, slit_width=ft_in.slit_width, n_slits=ft_in.n_slits,
                N_seg=ft_in.N_seg,
                circuit_mode=ft_in.circuit_mode,
                circuits=ft_in.circuits or [],
            )

        # Build MCHX spec
        mchx = None
        if req.hx_type == "MCHX":
            m_in = req.mchx_spec or MCHXSpecInput()
            mchx = MCHXSpec(
                W=m_in.W, H=m_in.H, D=m_in.D,
                Nr=m_in.Nr,
                ch_width=m_in.ch_width, ch_height=m_in.ch_height,
                ch_wall=m_in.ch_wall, n_ports=m_in.n_ports,
                tube_height=m_in.tube_height, tube_pitch=m_in.tube_pitch,
                louver_pitch=m_in.louver_pitch, louver_angle=m_in.louver_angle,
                fin_pitch=m_in.fin_pitch, fin_thickness=m_in.fin_thickness,
                fin_height=m_in.fin_height, k_fin=m_in.k_fin,
                N_seg=m_in.N_seg, Nt=m_in.Nt,
                passes=m_in.passes or [],
                pass_slabs=m_in.pass_slabs or [],
            )

        # --- Resolve V_air ↔ CMM ---
        # Frontal area: H × W for both FT and MCHX
        if req.hx_type == "FT":
            face_H = (req.ft_spec or FTSpecInput()).H
            face_W = (req.ft_spec or FTSpecInput()).W
        else:
            face_H = (req.mchx_spec or MCHXSpecInput()).H
            face_W = (req.mchx_spec or MCHXSpecInput()).W
        A_face = face_H * face_W  # [m²]

        if req.air_flow_mode == "CMM" and req.CMM is not None:
            # CMM [m³/min] → V_air [m/s]
            V_air_resolved = (req.CMM / 60.0) / A_face if A_face > 0 else 2.0
            CMM_resolved = req.CMM
        else:
            # V_air [m/s] → CMM [m³/min]
            V_air_resolved = req.V_air if req.V_air is not None else 2.0
            CMM_resolved = V_air_resolved * A_face * 60.0

        # Determine T_ref_in for single-phase entry
        T_ref_in_K = None
        if not req.two_phase_inlet and req.T_ref_in_C is not None:
            T_ref_in_K = req.T_ref_in_C + 273.15

        sim_input = SimulationInput(
            hx_type=req.hx_type,
            mode=req.mode,
            T_air_in=T_air_K,
            RH_in=req.RH_in,
            V_air=V_air_resolved,
            fluid=req.fluid,
            T_sat=T_sat_K,
            m_ref=req.m_ref,
            x_in=x_in,
            T_ref_in=T_ref_in_K,
            flow_arrangement=req.flow_arrangement,
            max_outer=req.max_outer,
            outer_tol_pct=req.outer_tol_pct,
            cf_j=req.cf_j,
            cf_f=req.cf_f,
            cf_hi=req.cf_hi,
            cf_dp_ref=req.cf_dp_ref,
            ft_spec=ft,
            mchx_spec=mchx,
        )

        # --- Correlation recommendation ---
        from models.geometry import FinTubeGeo, MCHXGeo
        from models.properties import MoistAirProperties
        air_props = MoistAirProperties()
        mu_air = air_props.mu_air(T_air_K)

        if req.hx_type == "FT":
            ft_resolved = ft or FinTubeSpec()
            geo_temp = FinTubeGeo.from_spec(ft_resolved)
            Dc = geo_temp.Dc
            rho_air = air_props.rho_air(T_air_K, 0.01)
            G_air = rho_air * V_air_resolved / geo_temp.sigma if geo_temp.sigma > 0 else 5.0
            Re_Dc_val = G_air * Dc / mu_air
            fin_type = ft_resolved.fin_type
            Nr = ft_resolved.Nr
            spec_vals = build_spec_values(ft_resolved, geo_temp, Re_Dc_val)
        else:
            mchx_resolved = mchx or MCHXSpec()
            geo_temp = MCHXGeo.from_spec(mchx_resolved)
            rho_air = air_props.rho_air(T_air_K, 0.01)
            G_air_mchx = rho_air * V_air_resolved / geo_temp.sigma if geo_temp.sigma > 0 else 5.0
            Re_Dc_val = G_air_mchx * mchx_resolved.louver_pitch / mu_air  # Re_Lp
            fin_type = "mchx"
            Nr = mchx_resolved.Nr
            spec_vals = build_mchx_spec_values(mchx_resolved, geo_temp, Re_Dc_val)

        rec = recommend_correlation(fin_type, Re_Dc_val, Nr, req.hx_type, spec_vals)

        # f-factor recommendation
        f_fin_type = fin_type if req.hx_type == "FT" else "mchx"
        f_rec = recommend_f_correlation(f_fin_type, Re_Dc_val, req.hx_type)

        # Apply user-selected or recommended correlation
        solver = HXSolver(sim_input)

        # j-factor selection
        if req.air_j_corr:
            solver.corr["air_j"] = req.air_j_corr
            user_val = validate_correlation(req.air_j_corr, spec_vals)
            rec["user_selected"] = req.air_j_corr
            rec["user_validation"] = user_val
        else:
            solver.corr["air_j"] = rec["recommended"]

        # f-factor selection
        if req.air_f_corr:
            solver._f_corr_id = req.air_f_corr
            f_rec["user_selected"] = req.air_f_corr
        else:
            solver._f_corr_id = f_rec["recommended"]

        # Refrigerant-side correlation selection
        if req.evap_corr:
            solver.corr["evap"] = req.evap_corr
        if req.cond_corr:
            solver.corr["cond"] = req.cond_corr
        if req.dp_ref_corr:
            solver.corr["dp_ref"] = req.dp_ref_corr

        result = solver.solve()

        # Build segment output
        seg_out = [
            SegmentOut(
                row=s.row, tube=s.tube, seg=s.seg,
                Q=round(s.Q, 2),
                Q_sensible=round(s.Q_sensible, 2),
                Q_latent=round(s.Q_latent, 2),
                T_wall=round(s.T_wall - 273.15, 2),
                T_air_local=round(s.T_air_local - 273.15, 2),
                T_ref=round(s.T_ref - 273.15, 2),
                x_ref=round(s.x_ref, 4),
                h_i=round(s.h_i, 1),
                h_o=round(s.h_o, 1),
                eta_o=round(s.eta_o, 4),
                dp_ref=round(s.dp_ref, 2),
                is_wet=s.is_wet,
                converged=s.converged,
            )
            for s in result.segments
        ]

        return SimResponse(
            Q_total=round(result.Q_total, 1),
            Q_sensible=round(result.Q_sensible, 1),
            Q_latent=round(result.Q_latent, 1),
            SHR=round(result.SHR, 4),
            T_air_out_C=round(result.T_air_out - 273.15, 2),
            W_air_out=round(result.W_air_out, 6),
            RH_out=round(result.RH_out, 4),
            x_ref_out=round(result.x_ref_out, 4),
            T_ref_out_C=round(result.T_ref_out - 273.15, 2),
            dp_air=round(result.dp_air, 1),
            dp_ref=round(result.dp_ref, 1),
            T_sat_C=round(T_sat_C_resolved, 2),
            P_sat_kPa=round(P_sat_kPa_resolved, 1),
            V_air=round(V_air_resolved, 3),
            CMM=round(CMM_resolved, 3),
            Re_Dc=round(Re_Dc_val, 0),
            row_Q=[round(q, 1) for q in result.row_Q],
            correlations_used=result.correlations_used,
            correlation_recommendation=rec,
            f_recommendation=f_rec,
            convergence=result.convergence,
            segments=seg_out,
            error=result.error,
        )

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Simulation error: {str(e)}\n{traceback.format_exc()}")


@app.get("/correlations")
def get_correlations(fin_type: str = "plain"):
    """Get available j-factor correlations and metadata for a fin type."""
    available = get_available_correlations(fin_type)
    details = {cid: AIRSIDE_CORRELATIONS[cid] for cid in available if cid in AIRSIDE_CORRELATIONS}
    return {"fin_type": fin_type, "available": available, "details": details}


@app.get("/f_correlations")
def get_f_correlations(fin_type: str = "plain", hx_type: str = "FT"):
    """Get available f-factor correlations and metadata for a fin type."""
    ft = fin_type if hx_type == "FT" else "mchx"
    available = get_available_f_correlations(ft)
    details = {cid: FSIDE_CORRELATIONS[cid] for cid in available if cid in FSIDE_CORRELATIONS}
    return {"fin_type": fin_type, "hx_type": hx_type, "available": available, "details": details}


@app.get("/circuit_presets")
def get_circuit_presets(Nr: int = 4, Nt: int = 4, flow: str = "counter"):
    """Get all available circuit preset patterns."""
    from models.geometry import generate_circuits
    presets = {}
    for mode in ["row_parallel", "serpentine_2", "serpentine_4", "single"]:
        try:
            circuits = generate_circuits(Nr, Nt, mode, flow)
            presets[mode] = circuits
        except:
            pass
    return {"Nr": Nr, "Nt": Nt, "flow": flow, "presets": presets}


@app.get("/ref_correlations")
def get_ref_correlations(mode: str = "evap", hx_type: str = "FT"):
    """Get available refrigerant-side correlations filtered by HX type."""
    available = get_available_ref_correlations(mode, hx_type)
    registry = REFSIDE_EVAP_CORRELATIONS if mode == "evap" else REFSIDE_COND_CORRELATIONS
    details = {cid: registry[cid] for cid in available}
    return {"mode": mode, "hx_type": hx_type, "available": available, "details": details}


@app.get("/dp_ref_correlations")
def get_dp_ref_correlations_ep(hx_type: str = "FT"):
    """Get available refrigerant-side pressure drop correlations."""
    available = get_available_dp_correlations(hx_type)
    details = {cid: REFSIDE_DP_CORRELATIONS[cid] for cid in available}
    return {"hx_type": hx_type, "available": available, "details": details}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
