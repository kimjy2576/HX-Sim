"""
Quick test — verify simulation engine runs correctly.
"""
import sys
sys.path.insert(0, ".")

from models.properties import RefrigerantProperties, MoistAirProperties
from models.geometry import FinTubeSpec, FinTubeGeo, MCHXSpec, MCHXGeo
from models.correlations import select_correlations, h_evap_chen1966, h_cond_shah1979
from models.solver import SimulationInput, HXSolver


def test_properties():
    print("=" * 50)
    print("1. Refrigerant Properties (R410A)")
    ref = RefrigerantProperties("R410A")
    P = ref.P_sat(280.15)  # 7°C
    print(f"   P_sat(7°C)   = {P/1000:.1f} kPa")
    print(f"   T_sat(P)     = {ref.T_sat(P)-273.15:.2f} °C")
    print(f"   h_fg         = {ref.h_fg(P)/1000:.1f} kJ/kg")
    print(f"   rho_l        = {ref.rho_l(P):.1f} kg/m³")
    print(f"   rho_v        = {ref.rho_v(P):.2f} kg/m³")
    print(f"   P_r          = {ref.P_r(P):.4f}")
    print(f"   Xtt(x=0.5)   = {ref.Xtt(0.5, P):.4f}")
    print("   ✅ PASS")


def test_geometry():
    print("=" * 50)
    print("2. FT-HX Geometry")
    spec = FinTubeSpec(Nr=4, Nt=12, Di=0.00822, Do=0.00952)
    geo = FinTubeGeo.from_spec(spec)
    print(f"   A_total = {geo.A_total:.4f} m²")
    print(f"   A_i     = {geo.A_i:.4f} m²")
    print(f"   σ       = {geo.sigma:.4f}")
    print(f"   Dh      = {geo.Dh*1000:.2f} mm")

    eta_fin, eta_o = geo.fin_efficiency_schmidt(100.0)
    print(f"   η_fin(h=100) = {eta_fin:.4f}")
    print(f"   η_o(h=100)   = {eta_o:.4f}")

    print("\n   MCHX Geometry")
    mspec = MCHXSpec()
    mgeo = MCHXGeo.from_spec(mspec)
    print(f"   Dh_ref  = {mgeo.Dh_ref*1000:.3f} mm")
    print(f"   A_i     = {mgeo.A_i:.4f} m²")
    print(f"   N_ch    = {mgeo.N_ch}")
    print("   ✅ PASS")


def test_correlations():
    print("=" * 50)
    print("3. Heat Transfer Correlations")
    ref = RefrigerantProperties("R410A")
    P = ref.P_sat(280.15)

    h_e = h_evap_chen1966(x=0.5, G=200, Di=0.00822, ref=ref, P=P)
    print(f"   Chen(1966) h_evap(x=0.5) = {h_e:.0f} W/m²K")

    h_c = h_cond_shah1979(x=0.5, G=200, Di=0.00822, ref=ref, P=P)
    print(f"   Shah(1979) h_cond(x=0.5) = {h_c:.0f} W/m²K")

    corr = select_correlations("FT", 0.00822, "plain", 0.0254, 0.022)
    print(f"   Auto-select(FT): {corr}")

    corr_m = select_correlations("MCHX", 0.0012, "louver")
    print(f"   Auto-select(MCHX): {corr_m}")
    print("   ✅ PASS")


def test_simulation():
    print("=" * 50)
    print("4. Full Simulation — FT Evaporator")
    inp = SimulationInput(
        hx_type="FT",
        mode="evap",
        T_air_in=308.15,  # 35°C
        RH_in=0.50,
        V_air=2.0,
        fluid="R410A",
        T_sat=280.15,     # 7°C
        m_ref=0.02,
        x_in=0.2,
        flow_arrangement="counter",
        ft_spec=FinTubeSpec(Nr=4, Nt=12, N_seg=5),
    )
    solver = HXSolver(inp)
    result = solver.solve()

    print(f"   Q_total   = {result.Q_total:.1f} W")
    print(f"   Q_sens    = {result.Q_sensible:.1f} W")
    print(f"   Q_lat     = {result.Q_latent:.1f} W")
    print(f"   SHR       = {result.SHR:.4f}")
    print(f"   T_air_out = {result.T_air_out-273.15:.2f} °C")
    print(f"   x_ref_out = {result.x_ref_out:.4f}")
    print(f"   ΔP_air    = {result.dp_air:.1f} Pa")
    print(f"   Segments  = {len(result.segments)}")
    print(f"   Row Q     = {[round(q,1) for q in result.row_Q]}")
    print(f"   Corr used = {result.correlations_used}")

    if result.error:
        print(f"   ⚠️ Error: {result.error}")
    else:
        print("   ✅ PASS")

    print("\n5. Full Simulation — MCHX Condenser")
    inp2 = SimulationInput(
        hx_type="MCHX",
        mode="cond",
        T_air_in=308.15,
        RH_in=0.40,
        V_air=2.5,
        fluid="R410A",
        T_sat=318.15,     # 45°C
        m_ref=0.03,
        x_in=0.95,
        flow_arrangement="counter",
        mchx_spec=MCHXSpec(N_tubes=40, n_slabs=1, N_seg=5),
    )
    solver2 = HXSolver(inp2)
    result2 = solver2.solve()

    print(f"   Q_total   = {result2.Q_total:.1f} W")
    print(f"   SHR       = {result2.SHR:.4f}")
    print(f"   T_air_out = {result2.T_air_out-273.15:.2f} °C")
    print(f"   x_ref_out = {result2.x_ref_out:.4f}")
    print(f"   Row Q     = {[round(q,1) for q in result2.row_Q]}")

    if result2.error:
        print(f"   ⚠️ Error: {result2.error}")
    else:
        print("   ✅ PASS")


if __name__ == "__main__":
    print("🔬 HX Simulator v7.9-L2 — Test Suite")
    print()
    test_properties()
    print()
    test_geometry()
    print()
    test_correlations()
    print()
    test_simulation()
    print()
    print("=" * 50)
    print("All tests complete!")
