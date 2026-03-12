import { useState, useCallback } from "react";
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, BarChart, Bar, Cell } from "recharts";

const API_URL = "http://localhost:8000";

// ============================================================
// Default specs
// ============================================================
const DEFAULT_FT = {
  W: 0.5, H: 0.3, D: 0.08,
  Do: 0.00952, Di: 0.00822,
  Pt: 0.0254, Pl: 0.022,
  Nr: 4, Nt: 12,
  layout: "staggered",
  FPI: 14, fin_thickness: 0.00012,
  fin_type: "plain", k_fin: 200,
  louver_pitch: 0.0017, louver_angle: 27,
  N_seg: 5,
};

const DEFAULT_MCHX = {
  W: 0.6, H: 0.4, D: 0.02,
  n_slabs: 1,
  ch_width: 0.001, ch_height: 0.0015,
  ch_wall: 0.0003, n_ports: 12,
  tube_height: 0.002, tube_pitch: 0.01,
  louver_pitch: 0.0013, louver_angle: 27,
  fin_pitch: 0.0014, fin_thickness: 0.0001,
  fin_height: 0.008, k_fin: 200,
  N_seg: 5, N_tubes: 40,
};

const REFRIGERANTS = ["R410A","R134a","R32","R290","R1234yf","R22","R407C","R404A","R513A","R454C"];

// ============================================================
// Styles
// ============================================================
const palette = {
  bg: "#0a0e17",
  surface: "#111827",
  card: "#1a2233",
  border: "#2a3a52",
  accent: "#38bdf8",
  accent2: "#818cf8",
  accent3: "#f472b6",
  warn: "#fbbf24",
  success: "#34d399",
  text: "#e2e8f0",
  textDim: "#94a3b8",
  textBright: "#f8fafc",
  evap: "#38bdf8",
  cond: "#f472b6",
};

// ============================================================
// Input Group Component
// ============================================================
function InputGroup({ label, children }) {
  return (
    <div style={{
      background: palette.card,
      border: `1px solid ${palette.border}`,
      borderRadius: 12,
      padding: "14px 16px",
      marginBottom: 12,
    }}>
      <div style={{
        fontSize: 11, fontWeight: 700, textTransform: "uppercase",
        letterSpacing: 1.5, color: palette.accent, marginBottom: 10,
      }}>{label}</div>
      {children}
    </div>
  );
}

function Field({ label, unit, value, onChange, type = "number", options, step, min, max }) {
  const st = {
    display: "flex", alignItems: "center", justifyContent: "space-between",
    marginBottom: 6, gap: 8,
  };
  const inputSt = {
    background: palette.bg, color: palette.text,
    border: `1px solid ${palette.border}`, borderRadius: 6,
    padding: "5px 8px", fontSize: 13, width: 100, textAlign: "right",
    outline: "none",
  };
  const selectSt = { ...inputSt, width: 120, textAlign: "left", cursor: "pointer" };

  if (options) {
    return (
      <div style={st}>
        <span style={{ fontSize: 12, color: palette.textDim }}>{label}</span>
        <select value={value} onChange={e => onChange(e.target.value)} style={selectSt}>
          {options.map(o => <option key={o} value={o}>{o}</option>)}
        </select>
      </div>
    );
  }

  return (
    <div style={st}>
      <span style={{ fontSize: 12, color: palette.textDim }}>
        {label} {unit && <span style={{ color: palette.accent, fontSize: 10 }}>[{unit}]</span>}
      </span>
      <input
        type={type} value={value} step={step || "any"} min={min} max={max}
        onChange={e => onChange(type === "number" ? parseFloat(e.target.value) || 0 : e.target.value)}
        style={inputSt}
      />
    </div>
  );
}

// ============================================================
// Results Cards
// ============================================================
function MetricCard({ label, value, unit, color = palette.accent }) {
  return (
    <div style={{
      background: palette.card, border: `1px solid ${palette.border}`,
      borderRadius: 10, padding: "12px 16px", flex: "1 1 140px", minWidth: 140,
      borderTop: `3px solid ${color}`,
    }}>
      <div style={{ fontSize: 10, color: palette.textDim, textTransform: "uppercase", letterSpacing: 1 }}>{label}</div>
      <div style={{ fontSize: 22, fontWeight: 700, color, marginTop: 4 }}>
        {value} <span style={{ fontSize: 11, color: palette.textDim }}>{unit}</span>
      </div>
    </div>
  );
}

// ============================================================
// Main App
// ============================================================
export default function HXSimulator() {
  // State
  const [hxType, setHxType] = useState("FT");
  const [mode, setMode] = useState("evap");
  const [airIn, setAirIn] = useState({ T: 35, RH: 0.5, V: 2.0 });
  const [refIn, setRefIn] = useState({ fluid: "R410A", T_sat: 7, m_ref: 0.02, x_in: 0.2 });
  const [flow, setFlow] = useState("counter");
  const [ftSpec, setFtSpec] = useState({ ...DEFAULT_FT });
  const [mchxSpec, setMchxSpec] = useState({ ...DEFAULT_MCHX });
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const [tab, setTab] = useState("overview");

  const updateAir = (k, v) => setAirIn(p => ({ ...p, [k]: v }));
  const updateRef = (k, v) => setRefIn(p => ({ ...p, [k]: v }));
  const updateFT = (k, v) => setFtSpec(p => ({ ...p, [k]: v }));
  const updateMCHX = (k, v) => setMchxSpec(p => ({ ...p, [k]: v }));

  // Run simulation
  const runSim = useCallback(async () => {
    setLoading(true);
    setError("");
    try {
      const body = {
        hx_type: hxType,
        mode,
        T_air_in_C: airIn.T,
        RH_in: airIn.RH,
        V_air: airIn.V,
        fluid: refIn.fluid,
        T_sat_C: refIn.T_sat,
        m_ref: refIn.m_ref,
        x_in: refIn.x_in,
        flow_arrangement: flow,
      };
      if (hxType === "FT") body.ft_spec = ftSpec;
      else body.mchx_spec = mchxSpec;

      const res = await fetch(`${API_URL}/simulate`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
      if (!res.ok) {
        const err = await res.json();
        throw new Error(err.detail || "Simulation failed");
      }
      const data = await res.json();
      setResult(data);
      setTab("overview");
    } catch (e) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  }, [hxType, mode, airIn, refIn, flow, ftSpec, mchxSpec]);

  // Prepare chart data
  const segChartData = result?.segments
    ?.filter(s => s.tube === 0)
    ?.map((s, i) => ({
      name: `R${s.row}-S${s.seg}`,
      Q: s.Q,
      T_wall: s.T_wall,
      h_i: s.h_i,
      x: s.x_ref,
      eta_o: s.eta_o,
    })) || [];

  const rowChartData = result?.row_Q?.map((q, i) => ({ name: `Row ${i}`, Q: q })) || [];

  const modeColor = mode === "evap" ? palette.evap : palette.cond;

  return (
    <div style={{
      background: palette.bg, color: palette.text, minHeight: "100vh",
      fontFamily: "'JetBrains Mono', 'SF Mono', 'Fira Code', monospace",
    }}>
      {/* Header */}
      <div style={{
        background: `linear-gradient(135deg, ${palette.surface} 0%, ${palette.card} 100%)`,
        borderBottom: `1px solid ${palette.border}`,
        padding: "16px 24px", display: "flex", alignItems: "center", gap: 16,
      }}>
        <div style={{ fontSize: 20, fontWeight: 800, letterSpacing: -0.5 }}>
          <span style={{ color: palette.accent }}>❄️</span> HX Simulator
          <span style={{ color: palette.textDim, fontSize: 11, marginLeft: 8 }}>v7.9-L2</span>
        </div>
        <div style={{ flex: 1 }} />
        <div style={{
          fontSize: 10, color: palette.textDim, textAlign: "right", lineHeight: 1.4,
        }}>
          Chen(1966) · Kim&Mudawar(2013/12) · Shah(1979)<br />
          Gnielinski(1976) · Wang(2000) · Chang&Wang(1997)
        </div>
      </div>

      <div style={{ display: "flex", minHeight: "calc(100vh - 60px)" }}>
        {/* ===== LEFT PANEL: Inputs ===== */}
        <div style={{
          width: 320, flexShrink: 0, padding: 16,
          borderRight: `1px solid ${palette.border}`,
          overflowY: "auto", maxHeight: "calc(100vh - 60px)",
        }}>
          {/* HX Type & Mode */}
          <InputGroup label="Configuration">
            <div style={{ display: "flex", gap: 6, marginBottom: 8 }}>
              {["FT", "MCHX"].map(t => (
                <button key={t} onClick={() => setHxType(t)} style={{
                  flex: 1, padding: "8px 0", borderRadius: 8, border: "none", cursor: "pointer",
                  fontWeight: 700, fontSize: 13, fontFamily: "inherit",
                  background: hxType === t ? palette.accent : palette.bg,
                  color: hxType === t ? palette.bg : palette.textDim,
                  transition: "all 0.2s",
                }}>{t === "FT" ? "🔧 Fin-Tube" : "⚡ MCHX"}</button>
              ))}
            </div>
            <div style={{ display: "flex", gap: 6 }}>
              {["evap", "cond"].map(m => (
                <button key={m} onClick={() => setMode(m)} style={{
                  flex: 1, padding: "8px 0", borderRadius: 8, border: "none", cursor: "pointer",
                  fontWeight: 700, fontSize: 13, fontFamily: "inherit",
                  background: mode === m ? modeColor : palette.bg,
                  color: mode === m ? palette.bg : palette.textDim,
                }}>{m === "evap" ? "❄️ Evaporator" : "🔥 Condenser"}</button>
              ))}
            </div>
            <div style={{ marginTop: 8 }}>
              <Field label="Flow" options={["counter", "parallel"]} value={flow} onChange={setFlow} />
            </div>
          </InputGroup>

          {/* Air Inlet */}
          <InputGroup label="🌀 Air Inlet">
            <Field label="T_air" unit="°C" value={airIn.T} onChange={v => updateAir("T", v)} />
            <Field label="RH" unit="-" value={airIn.RH} onChange={v => updateAir("RH", v)} step={0.05} min={0} max={1} />
            <Field label="V_face" unit="m/s" value={airIn.V} onChange={v => updateAir("V", v)} />
          </InputGroup>

          {/* Refrigerant */}
          <InputGroup label="🧊 Refrigerant">
            <Field label="Fluid" options={REFRIGERANTS} value={refIn.fluid} onChange={v => updateRef("fluid", v)} />
            <Field label="T_sat" unit="°C" value={refIn.T_sat} onChange={v => updateRef("T_sat", v)} />
            <Field label="ṁ_ref" unit="kg/s" value={refIn.m_ref} onChange={v => updateRef("m_ref", v)} step={0.001} />
            <Field label="x_in" unit="-" value={refIn.x_in} onChange={v => updateRef("x_in", v)} step={0.05} />
          </InputGroup>

          {/* Geometry */}
          {hxType === "FT" ? (
            <InputGroup label="📐 FT Geometry">
              <Field label="W" unit="m" value={ftSpec.W} onChange={v => updateFT("W", v)} />
              <Field label="H" unit="m" value={ftSpec.H} onChange={v => updateFT("H", v)} />
              <Field label="Do" unit="mm" value={Math.round(ftSpec.Do * 1000 * 100) / 100}
                     onChange={v => updateFT("Do", v / 1000)} />
              <Field label="Di" unit="mm" value={Math.round(ftSpec.Di * 1000 * 100) / 100}
                     onChange={v => updateFT("Di", v / 1000)} />
              <Field label="Pt" unit="mm" value={Math.round(ftSpec.Pt * 1000 * 10) / 10}
                     onChange={v => updateFT("Pt", v / 1000)} />
              <Field label="Pl" unit="mm" value={Math.round(ftSpec.Pl * 1000 * 10) / 10}
                     onChange={v => updateFT("Pl", v / 1000)} />
              <Field label="Nr" unit="rows" value={ftSpec.Nr} onChange={v => updateFT("Nr", parseInt(v) || 1)} />
              <Field label="Nt" unit="tubes" value={ftSpec.Nt} onChange={v => updateFT("Nt", parseInt(v) || 1)} />
              <Field label="FPI" unit="fins/in" value={ftSpec.FPI} onChange={v => updateFT("FPI", v)} />
              <Field label="Fin type" options={["plain","wavy","louver","slit"]} value={ftSpec.fin_type}
                     onChange={v => updateFT("fin_type", v)} />
              <Field label="N_seg" unit="per tube" value={ftSpec.N_seg} onChange={v => updateFT("N_seg", parseInt(v) || 1)} />
            </InputGroup>
          ) : (
            <InputGroup label="📐 MCHX Geometry">
              <Field label="W" unit="m" value={mchxSpec.W} onChange={v => updateMCHX("W", v)} />
              <Field label="H" unit="m" value={mchxSpec.H} onChange={v => updateMCHX("H", v)} />
              <Field label="ch_w" unit="mm" value={mchxSpec.ch_width * 1000}
                     onChange={v => updateMCHX("ch_width", v / 1000)} />
              <Field label="ch_h" unit="mm" value={mchxSpec.ch_height * 1000}
                     onChange={v => updateMCHX("ch_height", v / 1000)} />
              <Field label="n_ports" unit="-" value={mchxSpec.n_ports}
                     onChange={v => updateMCHX("n_ports", parseInt(v) || 1)} />
              <Field label="N_tubes" unit="-" value={mchxSpec.N_tubes}
                     onChange={v => updateMCHX("N_tubes", parseInt(v) || 1)} />
              <Field label="Lp" unit="mm" value={mchxSpec.louver_pitch * 1000}
                     onChange={v => updateMCHX("louver_pitch", v / 1000)} />
              <Field label="θ_louver" unit="°" value={mchxSpec.louver_angle}
                     onChange={v => updateMCHX("louver_angle", v)} />
              <Field label="Fp" unit="mm" value={mchxSpec.fin_pitch * 1000}
                     onChange={v => updateMCHX("fin_pitch", v / 1000)} />
              <Field label="n_slabs" unit="-" value={mchxSpec.n_slabs}
                     onChange={v => updateMCHX("n_slabs", parseInt(v) || 1)} />
              <Field label="N_seg" unit="per tube" value={mchxSpec.N_seg}
                     onChange={v => updateMCHX("N_seg", parseInt(v) || 1)} />
            </InputGroup>
          )}

          {/* Run Button */}
          <button onClick={runSim} disabled={loading} style={{
            width: "100%", padding: "14px 0", borderRadius: 10,
            border: "none", cursor: loading ? "wait" : "pointer",
            fontWeight: 800, fontSize: 15, fontFamily: "inherit",
            letterSpacing: 1,
            background: loading
              ? palette.textDim
              : `linear-gradient(135deg, ${palette.accent} 0%, ${palette.accent2} 100%)`,
            color: palette.bg,
            boxShadow: `0 4px 24px ${palette.accent}33`,
            transition: "all 0.3s",
          }}>
            {loading ? "⏳ Computing..." : "▶ RUN SIMULATION"}
          </button>

          {error && (
            <div style={{
              marginTop: 12, padding: 12, borderRadius: 8,
              background: "#7f1d1d33", border: "1px solid #ef4444",
              fontSize: 11, color: "#fca5a5", wordBreak: "break-all",
            }}>{error}</div>
          )}
        </div>

        {/* ===== RIGHT PANEL: Results ===== */}
        <div style={{ flex: 1, padding: 20, overflowY: "auto", maxHeight: "calc(100vh - 60px)" }}>
          {!result ? (
            <div style={{
              display: "flex", alignItems: "center", justifyContent: "center",
              height: "100%", flexDirection: "column", gap: 12,
            }}>
              <div style={{ fontSize: 48 }}>🔬</div>
              <div style={{ color: palette.textDim, fontSize: 14 }}>
                Configure parameters and click <b>RUN SIMULATION</b>
              </div>
              <div style={{ color: palette.textDim, fontSize: 11, maxWidth: 400, textAlign: "center", lineHeight: 1.6 }}>
                Level 2 Tube-Segment Model: Nr×Nt×N_seg 세그먼트별 독립 계산, T_wall 반복 수렴.<br />
                Backend requires FastAPI + CoolProp (localhost:8000).
              </div>
            </div>
          ) : (
            <>
              {/* Tabs */}
              <div style={{ display: "flex", gap: 4, marginBottom: 16 }}>
                {["overview", "segments", "charts"].map(t => (
                  <button key={t} onClick={() => setTab(t)} style={{
                    padding: "8px 16px", borderRadius: 8, border: "none",
                    cursor: "pointer", fontFamily: "inherit", fontSize: 12,
                    fontWeight: 700, textTransform: "uppercase", letterSpacing: 1,
                    background: tab === t ? modeColor : palette.card,
                    color: tab === t ? palette.bg : palette.textDim,
                  }}>{t}</button>
                ))}
              </div>

              {/* OVERVIEW TAB */}
              {tab === "overview" && (
                <>
                  <div style={{ display: "flex", flexWrap: "wrap", gap: 10, marginBottom: 16 }}>
                    <MetricCard label="Q_total" value={result.Q_total.toFixed(1)} unit="W" color={modeColor} />
                    <MetricCard label="Q_sensible" value={result.Q_sensible.toFixed(1)} unit="W" color={palette.accent} />
                    <MetricCard label="Q_latent" value={result.Q_latent.toFixed(1)} unit="W" color={palette.accent2} />
                    <MetricCard label="SHR" value={result.SHR.toFixed(3)} unit="-" color={palette.warn} />
                  </div>
                  <div style={{ display: "flex", flexWrap: "wrap", gap: 10, marginBottom: 16 }}>
                    <MetricCard label="T_air,out" value={result.T_air_out_C.toFixed(1)} unit="°C" color={palette.success} />
                    <MetricCard label="x_ref,out" value={result.x_ref_out.toFixed(3)} unit="-" color={palette.accent3} />
                    <MetricCard label="ΔP_air" value={result.dp_air.toFixed(1)} unit="Pa" color={palette.textDim} />
                    <MetricCard label="RH_out" value={(result.RH_out * 100).toFixed(1)} unit="%" color={palette.accent2} />
                  </div>

                  {/* Correlations Used */}
                  <div style={{
                    background: palette.card, border: `1px solid ${palette.border}`,
                    borderRadius: 10, padding: 14, marginBottom: 16,
                  }}>
                    <div style={{ fontSize: 11, color: palette.accent, fontWeight: 700, marginBottom: 8, letterSpacing: 1 }}>
                      CORRELATIONS USED
                    </div>
                    <div style={{ display: "flex", flexWrap: "wrap", gap: 8 }}>
                      {Object.entries(result.correlations_used).map(([k, v]) => (
                        <span key={k} style={{
                          background: palette.bg, padding: "4px 10px", borderRadius: 6,
                          fontSize: 11, color: palette.text,
                        }}>
                          <span style={{ color: palette.textDim }}>{k}:</span> {v}
                        </span>
                      ))}
                    </div>
                  </div>

                  {/* Row-by-Row Q */}
                  <div style={{
                    background: palette.card, border: `1px solid ${palette.border}`,
                    borderRadius: 10, padding: 14,
                  }}>
                    <div style={{ fontSize: 11, color: palette.accent, fontWeight: 700, marginBottom: 10, letterSpacing: 1 }}>
                      ROW Q DISTRIBUTION
                    </div>
                    <ResponsiveContainer width="100%" height={180}>
                      <BarChart data={rowChartData}>
                        <CartesianGrid strokeDasharray="3 3" stroke={palette.border} />
                        <XAxis dataKey="name" tick={{ fill: palette.textDim, fontSize: 11 }} />
                        <YAxis tick={{ fill: palette.textDim, fontSize: 11 }} />
                        <Tooltip
                          contentStyle={{ background: palette.surface, border: `1px solid ${palette.border}`, borderRadius: 8, fontSize: 12 }}
                          labelStyle={{ color: palette.text }}
                        />
                        <Bar dataKey="Q" radius={[4, 4, 0, 0]}>
                          {rowChartData.map((_, i) => (
                            <Cell key={i} fill={i === 0 ? modeColor : `${modeColor}${Math.max(30, 99 - i * 20).toString()}`} />
                          ))}
                        </Bar>
                      </BarChart>
                    </ResponsiveContainer>
                  </div>
                </>
              )}

              {/* SEGMENTS TAB */}
              {tab === "segments" && (
                <div style={{
                  background: palette.card, border: `1px solid ${palette.border}`,
                  borderRadius: 10, padding: 14, overflowX: "auto",
                }}>
                  <div style={{ fontSize: 11, color: palette.accent, fontWeight: 700, marginBottom: 10, letterSpacing: 1 }}>
                    SEGMENT DETAILS (Tube 0)
                  </div>
                  <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 11 }}>
                    <thead>
                      <tr style={{ borderBottom: `1px solid ${palette.border}` }}>
                        {["Row","Seg","Q [W]","T_wall [°C]","T_air [°C]","x_ref","h_i","h_o","η_o","Wet","Conv."].map(h => (
                          <th key={h} style={{ padding: "6px 8px", textAlign: "right", color: palette.textDim, fontWeight: 600 }}>
                            {h}
                          </th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {result.segments.filter(s => s.tube === 0).map((s, i) => (
                        <tr key={i} style={{
                          borderBottom: `1px solid ${palette.border}22`,
                          background: i % 2 === 0 ? "transparent" : `${palette.bg}55`,
                        }}>
                          <td style={{ padding: "5px 8px", textAlign: "right" }}>{s.row}</td>
                          <td style={{ padding: "5px 8px", textAlign: "right" }}>{s.seg}</td>
                          <td style={{ padding: "5px 8px", textAlign: "right", color: modeColor, fontWeight: 600 }}>{s.Q.toFixed(1)}</td>
                          <td style={{ padding: "5px 8px", textAlign: "right" }}>{s.T_wall.toFixed(1)}</td>
                          <td style={{ padding: "5px 8px", textAlign: "right" }}>{s.T_air_local.toFixed(1)}</td>
                          <td style={{ padding: "5px 8px", textAlign: "right" }}>{s.x_ref.toFixed(3)}</td>
                          <td style={{ padding: "5px 8px", textAlign: "right" }}>{s.h_i.toFixed(0)}</td>
                          <td style={{ padding: "5px 8px", textAlign: "right" }}>{s.h_o.toFixed(0)}</td>
                          <td style={{ padding: "5px 8px", textAlign: "right" }}>{s.eta_o.toFixed(3)}</td>
                          <td style={{ padding: "5px 8px", textAlign: "center" }}>{s.is_wet ? "💧" : "—"}</td>
                          <td style={{ padding: "5px 8px", textAlign: "center", color: s.converged ? palette.success : "#ef4444" }}>
                            {s.converged ? "✓" : "✗"}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}

              {/* CHARTS TAB */}
              {tab === "charts" && (
                <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
                  {/* Q distribution */}
                  <div style={{ background: palette.card, border: `1px solid ${palette.border}`, borderRadius: 10, padding: 14 }}>
                    <div style={{ fontSize: 11, color: palette.accent, fontWeight: 700, marginBottom: 10, letterSpacing: 1 }}>
                      Q PER SEGMENT (Tube 0)
                    </div>
                    <ResponsiveContainer width="100%" height={220}>
                      <BarChart data={segChartData}>
                        <CartesianGrid strokeDasharray="3 3" stroke={palette.border} />
                        <XAxis dataKey="name" tick={{ fill: palette.textDim, fontSize: 10 }} />
                        <YAxis tick={{ fill: palette.textDim, fontSize: 10 }} />
                        <Tooltip contentStyle={{ background: palette.surface, border: `1px solid ${palette.border}`, borderRadius: 8, fontSize: 11 }} />
                        <Bar dataKey="Q" fill={modeColor} radius={[3, 3, 0, 0]} />
                      </BarChart>
                    </ResponsiveContainer>
                  </div>

                  {/* T_wall & x profile */}
                  <div style={{ background: palette.card, border: `1px solid ${palette.border}`, borderRadius: 10, padding: 14 }}>
                    <div style={{ fontSize: 11, color: palette.accent, fontWeight: 700, marginBottom: 10, letterSpacing: 1 }}>
                      T_WALL & QUALITY PROFILE
                    </div>
                    <ResponsiveContainer width="100%" height={220}>
                      <LineChart data={segChartData}>
                        <CartesianGrid strokeDasharray="3 3" stroke={palette.border} />
                        <XAxis dataKey="name" tick={{ fill: palette.textDim, fontSize: 10 }} />
                        <YAxis yAxisId="T" tick={{ fill: palette.textDim, fontSize: 10 }} />
                        <YAxis yAxisId="x" orientation="right" tick={{ fill: palette.textDim, fontSize: 10 }} domain={[0, 1.2]} />
                        <Tooltip contentStyle={{ background: palette.surface, border: `1px solid ${palette.border}`, borderRadius: 8, fontSize: 11 }} />
                        <Legend wrapperStyle={{ fontSize: 11 }} />
                        <Line yAxisId="T" type="monotone" dataKey="T_wall" stroke={palette.warn} strokeWidth={2} dot={{ r: 3 }} name="T_wall [°C]" />
                        <Line yAxisId="x" type="monotone" dataKey="x" stroke={palette.accent3} strokeWidth={2} dot={{ r: 3 }} name="x_ref [-]" />
                      </LineChart>
                    </ResponsiveContainer>
                  </div>

                  {/* h_i profile */}
                  <div style={{ background: palette.card, border: `1px solid ${palette.border}`, borderRadius: 10, padding: 14 }}>
                    <div style={{ fontSize: 11, color: palette.accent, fontWeight: 700, marginBottom: 10, letterSpacing: 1 }}>
                      h_i & η_o PROFILE
                    </div>
                    <ResponsiveContainer width="100%" height={220}>
                      <LineChart data={segChartData}>
                        <CartesianGrid strokeDasharray="3 3" stroke={palette.border} />
                        <XAxis dataKey="name" tick={{ fill: palette.textDim, fontSize: 10 }} />
                        <YAxis yAxisId="h" tick={{ fill: palette.textDim, fontSize: 10 }} />
                        <YAxis yAxisId="eta" orientation="right" tick={{ fill: palette.textDim, fontSize: 10 }} domain={[0, 1]} />
                        <Tooltip contentStyle={{ background: palette.surface, border: `1px solid ${palette.border}`, borderRadius: 8, fontSize: 11 }} />
                        <Legend wrapperStyle={{ fontSize: 11 }} />
                        <Line yAxisId="h" type="monotone" dataKey="h_i" stroke={palette.accent} strokeWidth={2} dot={{ r: 3 }} name="h_i [W/m²K]" />
                        <Line yAxisId="eta" type="monotone" dataKey="eta_o" stroke={palette.success} strokeWidth={2} dot={{ r: 3 }} name="η_o [-]" />
                      </LineChart>
                    </ResponsiveContainer>
                  </div>
                </div>
              )}
            </>
          )}
        </div>
      </div>
    </div>
  );
}
