# ❄️ HX Simulator v7.9-L2

**열교환기 해석 시뮬레이터 — FT-HX / MCHX Level 2 Tube-Segment Model**

> Nr × Nt × N_seg 세그먼트별 독립 계산, T_wall 반복 수렴  
> Chen(1966) · Kim&Mudawar(2013/12) · Shah(1979) · Gnielinski(1976) · Wang(2000) · Chang&Wang(1997)

---

## 📐 Architecture

```
hx-simulator/
├── backend/                 # FastAPI + CoolProp
│   ├── models/
│   │   ├── properties.py    # CoolProp 물성 래퍼 (냉매 + 습공기)
│   │   ├── geometry.py      # FT-HX / MCHX 기하 모델
│   │   ├── correlations.py  # 열전달 상관식 전체
│   │   └── solver.py        # T_wall 반복 수렴 솔버
│   ├── main.py              # FastAPI 엔드포인트
│   ├── Dockerfile
│   └── requirements.txt
├── frontend/
│   └── HXSimulator.jsx      # React UI (Recharts 시각화)
├── docker-compose.yml
└── README.md
```

## 🔧 Features

| 구분 | 지원 항목 |
|------|-----------|
| **HX 타입** | Fin-Tube (FT) / Micro-Channel (MCHX) |
| **모드** | 증발기 / 응축기 |
| **공기측** | Wang(2000) plain/wavy/louver/slit, Chang&Wang(1997) |
| **냉매측 (증발)** | Chen(1966) FT, Kim&Mudawar(2013) MCHX |
| **냉매측 (응축)** | Shah(1979) FT, Kim&Mudawar(2012) MCHX |
| **단상** | Gnielinski(1976) — 과열증기 + 과냉액체 |
| **전이 블렌딩** | x=0.90~1.05 (증발→과열), x=-0.05~0 (과냉→이상) |
| **습면 모델** | b@T_fin_avg 반복수렴, Schmidt 등가원형핀 |
| **핀효율** | Schmidt(FT), 직선핀(MCHX), 습면 보정 |
| **유동배열** | Counter / Parallel |
| **냉매** | R410A, R134a, R32, R290, R1234yf, R22, R407C 등 |

## 🚀 Quick Start

### 1. Backend (로컬)

```bash
cd backend
pip install -r requirements.txt
python main.py
# → http://localhost:8000/docs (Swagger UI)
```

### 2. Docker

```bash
docker-compose up --build
# → http://localhost:8000
```

### 3. Frontend (Claude.ai Artifact)

`frontend/HXSimulator.jsx`를 Claude.ai에서 React Artifact로 직접 실행 가능.  
백엔드(localhost:8000)가 실행 중이어야 동작함.

---

## 🌐 배포 (Render.com)

### Backend 배포

1. GitHub 레포에 push
2. [Render.com](https://render.com) → New Web Service
3. 설정:
   - **Build Command**: `pip install -r backend/requirements.txt`
   - **Start Command**: `cd backend && uvicorn main:app --host 0.0.0.0 --port $PORT`
   - **Docker**: 또는 `backend/Dockerfile` 사용

4. 배포 완료 후 URL 획득 (예: `https://hx-sim.onrender.com`)
5. `HXSimulator.jsx`의 `API_URL`을 해당 URL로 변경

---

## 📊 API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| GET | `/` | 서버 상태 |
| GET | `/refrigerants` | 지원 냉매 목록 |
| POST | `/simulate` | 시뮬레이션 실행 |
| GET | `/docs` | Swagger UI |

### POST `/simulate` Example

```json
{
  "hx_type": "FT",
  "mode": "evap",
  "T_air_in_C": 35,
  "RH_in": 0.5,
  "V_air": 2.0,
  "fluid": "R410A",
  "T_sat_C": 7,
  "m_ref": 0.02,
  "x_in": 0.2,
  "flow_arrangement": "counter",
  "ft_spec": {
    "Nr": 4, "Nt": 12, "N_seg": 5,
    "Di": 0.00822, "Do": 0.00952,
    "Pt": 0.0254, "Pl": 0.022,
    "FPI": 14, "fin_type": "plain"
  }
}
```

---

## 🔬 상관식 자동 선택 엔진

```
select_correlations(hx_type, Di, fin_type, Pt, Pl)

조건              공기 j       증발 h_i        응축 h_i
─────────────────────────────────────────────────────
FT plain         Wang(2000)   Chen(1966)      Shah(1979)
FT wavy          Wang(1999)   Chen(1966)      Shah(1979)
FT louver        Wang(1999)   Chen(1966)      Shah(1979)
MCHX (D<3mm)     C&W(1997)    K&M(2013)       K&M(2012)
```

---

## ⚠️ 한계 및 향후 과제

- [ ] 냉매 사이클 연동 (현재: T_sat 고정)
- [ ] 팬 P-Q 연동 커브
- [ ] CoolProp 물성 캐싱 (속도 10× 개선)
- [ ] Chang&Wang(2006) 91샘플 일반화
- [ ] 습면 f-factor 보정

---

## 📖 References

1. Chen, J.C. (1966) "Correlation for boiling heat transfer..." I&EC Process Design
2. Kim, S.M. & Mudawar, I. (2013) IJHMT 58:718-734
3. Kim, S.M. & Mudawar, I. (2012) IJHMT 55:3246-3261
4. Shah, M.M. (1979) ASHRAE Trans. 85:202-211
5. Gnielinski, V. (1976) Int. Chemical Eng. 16:359-368
6. Wang, C.C. et al. (2000) IJHMT 43(15):2693-2700
7. Chang, Y.J. & Wang, C.C. (1997) IJHMT 40(3):533-544

---

**License**: MIT
