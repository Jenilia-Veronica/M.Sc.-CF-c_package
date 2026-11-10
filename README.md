# GBM Monte Carlo Observatory

**Stock price simulation using Geometric Brownian Motion: a Monte Carlo engine written in C, with an interactive web dashboard.**

The C program does *all* of the mathematics and doubles as a tiny local web server. The browser page only draws what the C engine sends back as JSON.

> M.Sc. Computational Finance, AMCS Department,
> Author: Jenilia Veronica J

---

## Features

- **Exact GBM path simulation** using the closed-form solution from Itô's lemma, so there is no time-discretisation error
- **Fan chart** of possible futures: 5/25/50/75/95 percentile bands plus 60 sample paths
- **Terminal price distribution** compared with the theoretical log-normal density
- **Risk measures**: 95% Value-at-Risk (VaR) and Expected Shortfall (ES)
- **European option pricing** (call and put) by antithetic-variate Monte Carlo, validated against the Black-Scholes formula and put-call parity
- **Convergence plot** showing the Monte Carlo price approaching Black-Scholes with error shrinking like 1/√N
- **Live controls** for drift, volatility, horizon, number of paths and strike, with animated charts

---

## Mathematical model

Stock prices follow the stochastic differential equation

```
dS = μ S dt + σ S dW
```

where μ is the expected return, σ is the volatility and W is a standard Brownian motion. The exact solution used by the engine is

```
S(t + Δt) = S(t) · exp( (μ − ½σ²)Δt + σ√Δt · Z ),   Z ~ N(0, 1)
```

Option prices are computed under the risk-neutral measure (μ replaced by r) and compared with Black-Scholes:

```
C = S0·N(d1) − K·e^(−rT)·N(d2)
P = K·e^(−rT)·N(−d2) − S0·N(−d1)
```

---

## Implementation details

| Component | Method |
|---|---|
| Uniform random numbers | xoshiro256** (seeded via SplitMix64) |
| Normal random numbers | Box-Muller transform |
| Price paths | Exact GBM, 120 time steps per path |
| Percentile bands | Sorting the simulated prices at each time step |
| VaR / ES | 5th-percentile terminal price and mean of the worst 5% |
| Option pricing | 30,000 antithetic pairs (Z and −Z) |
| Convergence study | Up to 100,000 plain Monte Carlo draws, logged at geometric intervals |
| Web server | Raw sockets (Winsock on Windows, POSIX sockets elsewhere), no external libraries |

**Fixed constants:** S0 = 100, risk-free rate r = 5%.

---

## Repository structure

```
M.Sc.-CF-c_package/
├── README.md
└── package2/
    ├── gbm_server.c      # C engine: GBM, risk metrics, options, HTTP server
    ├── index.html        # Dashboard: draws charts from the engine's JSON
    ├── gbmserver.exe     # Pre-built Windows executable
    └── c-abstract.pdf    # One-page project abstract
```

---

## Getting started

`index.html` must be in the same folder as the executable, and the program must be run from that folder.

### Windows (quickest)

Open the `package2` folder and double-click `gbmserver.exe`, or run it from a terminal:

```
cd package2
gbmserver.exe
```

### Build from source

**Linux / macOS**

```
cd package2
gcc -O2 -o gbmserver gbm_server.c -lm
./gbmserver
```

**Windows (MinGW)**

```
cd package2
gcc -O2 -o gbmserver.exe gbm_server.c -lm -lws2_32
gbmserver.exe
```

The program starts a server on `http://localhost:8080` (it tries ports 8080 to 8099 if 8080 is busy) and attempts to open your browser automatically. If it doesn't, open the address printed in the terminal. Press `Ctrl+C` to stop.

The server listens on the loopback interface only, so it is reachable from your own computer and nowhere else.

---

## Using the dashboard

| Control | Meaning | Dashboard range |
|---|---|---|
| Expected return μ | Annual drift | −10% to 40% |
| Volatility σ | Annual volatility | 5% to 80% |
| Horizon T | Years simulated | 0.25 to 3 |
| Paths N | Number of simulated price paths | 200 to 10,000 |
| Option strike K | Strike for the call and put | 60 to 150 |

Moving any slider re-runs the C engine instantly. The **Replay with a new random future** button draws a fresh random seed.

The cards show simulated results next to their theoretical values (mean, standard deviation, probability of profit), plus VaR, Expected Shortfall, and the Monte Carlo versus Black-Scholes option prices.

---

## API

The engine exposes a single JSON endpoint, which you can also call directly:

```
GET /api/sim?mu=0.10&sg=0.20&T=1&N=3000&K=105&seed=1
```

| Parameter | Default | Allowed range |
|---|---|---|
| `mu` | 0.10 | −0.5 to 1.0 |
| `sg` | 0.20 | 0.01 to 2.0 |
| `T` | 1.0 | 0.1 to 10 |
| `N` | 3000 | 100 to 20,000 |
| `K` | 105 | 1 to 1000 |
| `seed` | 1 | any non-negative integer |

Out-of-range values are clamped. The response contains `bands`, `paths`, `hist`, `st` (summary statistics), `opt` (option prices), `conv` (convergence data) and `ms` (compute time). The same seed always reproduces the same result.

`GET /` serves `index.html`.

---

## Requirements

- A C compiler (GCC or MinGW) with the math library, only if building from source
- A modern web browser

No third-party libraries are needed.

---

## Author

**Jenilia Veronica J**,
M.Sc. Computational Finance, AMCS Department, PSG TECH.
