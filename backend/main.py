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
    louver_pitch: float = 0.0017
    louver_angle: float = 27.0
    N_seg: int = 5


class MCHXSpecInput(BaseModel):
    W: float = 0.6
    H: float = 0.4
    D: float = 0.020
    n_slabs: int = 1
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
    N_tubes: int = 40


class SimRequest(BaseModel):
    hx_type: Literal["FT", "MCHX"] = "FT"
    mode: Literal["evap", "cond"] = "evap"

    # Air
    T_air_in_C: float = Field(35.0, description="Air inlet temperature [°C]")
    RH_in: float = Field(0.50, ge=0, le=1, description="Relative humidity [-]")
    V_air: float = Field(2.0, description="Face velocity [m/s]")

    # Refrigerant
    fluid: str = "R410A"
    T_sat_C: float = Field(7.0, description="Saturation temperature [°C]")
    m_ref: float = Field(0.02, description="Mass flow rate [kg/s]")
    x_in: float = Field(0.2, description="Inlet quality")

    # Flow
    flow_arrangement: Literal["counter", "parallel"] = "counter"

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
    row_Q: List[float]
    correlations_used: Dict
    segments: List[SegmentOut]
    error: str = ""


# ============================================================
# Endpoints
# ============================================================

STATIC_DIR = Path(__file__).parent / "static"


@app.get("/")
def root():
    """Serve frontend UI."""
    return FileResponse(STATIC_DIR / "index.html")


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
        T_sat_K = req.T_sat_C + 273.15

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
                louver_pitch=ft_in.louver_pitch, louver_angle=ft_in.louver_angle,
                N_seg=ft_in.N_seg,
            )

        # Build MCHX spec
        mchx = None
        if req.hx_type == "MCHX":
            m_in = req.mchx_spec or MCHXSpecInput()
            mchx = MCHXSpec(
                W=m_in.W, H=m_in.H, D=m_in.D,
                n_slabs=m_in.n_slabs,
                ch_width=m_in.ch_width, ch_height=m_in.ch_height,
                ch_wall=m_in.ch_wall, n_ports=m_in.n_ports,
                tube_height=m_in.tube_height, tube_pitch=m_in.tube_pitch,
                louver_pitch=m_in.louver_pitch, louver_angle=m_in.louver_angle,
                fin_pitch=m_in.fin_pitch, fin_thickness=m_in.fin_thickness,
                fin_height=m_in.fin_height, k_fin=m_in.k_fin,
                N_seg=m_in.N_seg, N_tubes=m_in.N_tubes,
            )

        sim_input = SimulationInput(
            hx_type=req.hx_type,
            mode=req.mode,
            T_air_in=T_air_K,
            RH_in=req.RH_in,
            V_air=req.V_air,
            fluid=req.fluid,
            T_sat=T_sat_K,
            m_ref=req.m_ref,
            x_in=req.x_in,
            flow_arrangement=req.flow_arrangement,
            ft_spec=ft,
            mchx_spec=mchx,
        )

        solver = HXSolver(sim_input)
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
            row_Q=[round(q, 1) for q in result.row_Q],
            correlations_used=result.correlations_used,
            segments=seg_out,
            error=result.error,
        )

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Simulation error: {str(e)}\n{traceback.format_exc()}")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
