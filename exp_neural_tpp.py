"""ARM C (exploratory, first run): a RECAST-style neural temporal point process vs temporal ETAS.

Question
--------
Does a flexible autoregressive neural TPP (Dascher-Cousineau, Shchur, Brodsky & Guennemann
2023, GRL, "Using deep learning for flexible and scalable earthquake forecasting" -- RECAST)
carry information about *when the next earthquake happens* that temporal ETAS does not?
The only admissible answer is incremental bits/event over a frozen ETAS baseline on a
temporal hold-back, per HYPOTHESIS_LEDGER S-3.

Everything below is FROZEN before the first run.

Split (frozen)
--------------
  full span      = [t_first, t_last] of the catalogue at the magnitude floor
  EXPLORATION    = first 70% of the span                (engine/splits.py rule; the last 30%
                                                         is HOLDOUT and is never read here)
  within exploration, by TIME fraction of the exploration span:
      train      = first 60%
      validation = next 10%
      test       = last 30%
  Asserted in code: min(test timestamp) > max(train timestamp), and the encoder is strictly
  causal (state after consuming event i depends only on events <= i; the density for tau_{i+1}
  is emitted from that state).

Models (all time-only; magnitudes are INPUTS, never targets)
------------------------------------------------------------
  1. POISSON        lambda = constant = train-window event rate.
  2. ETAS           lambda(t) = mu + sum_{t_i<t} K 10^{alpha (M_i - M0)} (t - t_i + c)^{-p}
                    -- the exp_h_etas.py form, identical symbols, t in DAYS.
                    Fitted by MLE on the train window (L-BFGS-B on log-parameters, the same
                    multi-start seeds and bounds as exp_h_etas.py). The value and gradient are
                    computed with a torch port that is asserted (in this file and in
                    tests/test_neural_tpp.py) to agree with exp_h_etas.etas_ll to <1e-6
                    relative on a subsample; the torch port exists only because the numpy
                    version costs ~36 s per gradient evaluation at QTM scale.
  3. ETAS+DIURNAL   lambda(t) = g(hour(t)) * lambda_ETAS(t), g piecewise-constant on 24
                    local-solar-hour bins (the exact-integral equivalent of a Fourier series
                    truncated at the hourly Nyquist), fitted on train by its exact MLE
                    g_b = n_b / I_b with I_b the ETAS integral over the train times in bin b.
                    Local solar hour uses ONE reference longitude (the catalogue median), so
                    that lambda is a well-defined function of t alone.
  4. NEURAL TPP     GRU encoder over the event history; per-event inputs are
                    [z(log tau_i), z(M_i - M0), sin(2 pi h_i / 24), cos(2 pi h_i / 24)]
                    (z = train-window standardisation; the time-of-day pair is a DECLARED
                    nuisance input and is ablated). The state after event i emits a
                    K-component log-normal mixture for tau_{i+1} = t_{i+1} - t_i.
                    Trained with Adam, early stopping on validation LL.

Scoring (frozen) -- the SAME quantity for every model
-----------------------------------------------------
For test events i = te_lo .. N-1 (t_i in the test window), each model contributes
    log f(tau_i | history) = log lambda(t_i) - Integral_{t_{i-1}}^{t_i} lambda(u) du
which is a proper per-event inter-event-time log-density in units of 1/day. The intervals
[t_{i-1}, t_i] tile [t_{te_lo - 1}, t_{N-1}] exactly, so the total is a proper TPP
log-likelihood over the SAME events and the SAME window for every model, and -- the useful
part -- it decomposes EXACTLY per event, so every subset statistic below is a plain mean and
needs no apportionment of a global integral term.

    bits/event of A over B  =  mean_i [ log f_A(tau_i) - log f_B(tau_i) ] / ln 2

ETAS is scored WALK-FORWARD exactly as exp_h_etas.py: parameters frozen from the train fit,
history = all prior events (train + validation + earlier test).

PRIMARY statistic (frozen): bits/event of the neural TPP over ETAS+DIURNAL on the QTM test
window. Success rule (frozen): > 0.05 bits/event. The gain over raw ETAS is reported
alongside, never instead. Also reported: the gain restricted to test events with M >= 2.5
(where the day/night detection artifact is absent), and the gain over Poisson.

CALIBRATION CONTROL (mandatory, run before anything is believed)
----------------------------------------------------------------
An ETAS catalogue of the same size is simulated from the fitted train parameters (branching
process, GR magnitudes at the train b-value) and pushed through the IDENTICAL pipeline. On
ETAS-generated data the neural TPP must not beat ETAS by more than ~0.01 bits/event. If it
does, there is leakage or a scoring mismatch, and no real-data number is reported.

ARTIFACTS THAT COULD FAKE A WIN (PLAYBOOK standing rule 7 -- named before the run)
----------------------------------------------------------------------------------
1. FUTURE LEAKAGE. The catalogue is one long sequence; any non-causal encoder, any
   normalisation statistic computed off the training window, or any evaluation pass that
   lets a chunk see its own future would manufacture arbitrary bits. Guards: strict
   causality assertion on the encoder; all standardisation constants from train only; the
   evaluation pass is a single left-to-right run with the hidden state carried forward.
2. THE DAY/NIGHT DETECTION ARTIFACT. ~2% of sub-M2.5 events are displaced into local night
   by detection, not by the Earth. A flexible model WILL learn it. That is observer skill,
   not Earth skill. This is why the headline baseline is ETAS+DIURNAL (which is handed the
   same artifact) and why the M >= 2.5 restriction is reported.
3. COMPLETENESS DRIFT within 2008-2017. QTM is template-matched and far more uniform than a
   network catalogue, but its Mc is not constant in space or time; a model that learns
   "detection got better/worse" earns bits that are not physics.
4. ETAS MIS-SPECIFICATION BEING TRIVIALLY BEATEN. A single-kernel, spatially-blind,
   time-homogeneous ETAS is a weak straw man. Any win must be read as "ETAS as specified in
   exp_h_etas.py is missing something", not as "a new predictor of earthquakes".
5. TIED TIMESTAMPS. tau = 0 has infinite density under any continuous model, and a neural
   density with a numerical floor on tau would harvest unbounded "bits" from timestamp
   resolution rather than from physics. Events with tau == 0 are therefore DROPPED from
   scoring for EVERY model and the count is reported (1 in QTM at M>=1.0, 6 in SCSN at
   M>=2.5, out of 128k / 44k). Network inputs clamp tau at TAU_FLOOR_DAYS ~ 1 ms.
6. TRUNCATION MISMATCH. The Omori kernel is truncated at W days for tractability. The SAME W
   is used for every ETAS-family model; the truncation cost in bits/event is measured against
   an untruncated re-score and reported.

Attribution (section 7 of the brief) runs ONLY if the neural model wins on real data, and is
EXPLORATORY: it is not gated, not corrected for multiplicity, and claims nothing.

Outputs: results_neural_tpp.json, data/neural_tpp/*.pt
Run: python -u exp_neural_tpp.py
"""

from __future__ import annotations

import json
import math
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.optimize import minimize

import exp_h_etas as eh

HERE = Path(__file__).resolve().parent
OUT = HERE / "results_neural_tpp.json"
MODELDIR = HERE / "data" / "neural_tpp"
LN2 = math.log(2.0)
LN10 = math.log(10.0)

# ----------------------------------------------------------------- frozen configuration
SEED = 20260902

EXPLORE_FRAC = 0.70          # engine/splits.py
TRAIN_FRAC = 0.60            # of the exploration span
VAL_FRAC = 0.10
TEST_FRAC = 0.30

BURN_IN_DAYS = 365.0         # history-only lead-in before the fit/target window (as exp_h)
W_FIT_DAYS = 365.0           # Omori truncation for the train MLE
W_SCORE_DAYS = 1000.0        # Omori truncation for scoring (as exp_h's fit window); the
                             # untruncated correction is measured and reported
ALPHA_BOUNDS = eh.ALPHA_BOUNDS
P_BOUNDS = eh.P_BOUNDS
STARTS = eh.STARTS
MAXITER = 150

N_HOUR_BINS = 24
TAU_FLOOR_DAYS = 1e-8      # ~1 ms: the catalogue timestamp resolution

# neural TPP
NN = dict(
    hidden=128, layers=1, n_mix=32, dropout=0.0,
    lr=1e-3, batch=32, chunk_len=512, warmup_mask=128,
    grad_clip=1.0, max_steps=6000, eval_every=100, patience=12,
    max_train_minutes=35.0,
)
NN_CPU = dict(NN, hidden=64, n_mix=16, batch=16, chunk_len=256, warmup_mask=64,
              max_steps=1500, eval_every=100, patience=6, max_train_minutes=25.0)

GPU_BUDGET_BYTES = 6 * 1024**3   # hard cap on what we take from a shared card

SUCCESS_BITS = 0.05
CONTROL_TOL_BITS = 0.01

CATALOGS = {
    "QTM": dict(file="QTM_12dev.txt", M0=1.0, box=None),
    "SCSN": dict(file="SCSN_original_catalog.txt", M0=2.5,
                 box=dict(lat=[31.5, 38.0], lon=[-122.0, -113.5])),
}


# ----------------------------------------------------------------- device
def pick_device():
    """Choose a device under a hard memory budget. Records cudaMemGetInfo for every card."""
    info = {"cuda_available": bool(torch.cuda.is_available()), "devices": []}
    best, best_free = None, 0
    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            free, total = torch.cuda.mem_get_info(i)
            info["devices"].append({"index": i, "name": torch.cuda.get_device_name(i),
                                    "free_bytes": int(free), "total_bytes": int(total),
                                    "free_MiB": int(free // 2**20)})
            if free > best_free:
                best, best_free = i, free
    if best is not None and best_free >= GPU_BUDGET_BYTES:
        dev = torch.device(f"cuda:{best}")
        info["chosen"] = str(dev)
        info["reason"] = f"free {best_free // 2**20} MiB >= budget {GPU_BUDGET_BYTES // 2**20} MiB"
    else:
        dev = torch.device("cpu")
        info["chosen"] = "cpu"
        info["reason"] = (f"best free {best_free // 2**20} MiB < budget "
                          f"{GPU_BUDGET_BYTES // 2**20} MiB" if best is not None
                          else "no CUDA device")
    info["budget_bytes"] = GPU_BUDGET_BYTES
    return dev, info


# ----------------------------------------------------------------- ETAS in torch
class TorchETAS:
    """The exp_h_etas.py intensity, ported to torch. Same symbols, same formulas, t in days.

    Only two operations are needed:
      pair_S(query_times, W)        -> sum_{t_j < t} w_j (t - t_j + c)^{-p}   at query times
      seg_A(lo, hi, W)              -> Integral_{lo}^{hi} lambda(u) du        for many intervals
    Everything else (per-event log f, the diurnal factor, the fit) is built from these.
    """

    def __init__(self, t, m, device, dtype=torch.float64, budget=GPU_BUDGET_BYTES):
        self.device = device
        self.dtype = dtype
        self.t = torch.as_tensor(np.asarray(t, float), dtype=dtype, device=device)
        self.m = torch.as_tensor(np.asarray(m, float), dtype=dtype, device=device)
        self.n = self.t.numel()
        self.t_np = np.asarray(t, float)
        self.budget = budget

    def _chunk(self, n_src):
        """Rows per chunk so that ~8 temporaries of (rows x n_src) fit the budget."""
        itemsize = 8 if self.dtype == torch.float64 else 4
        rows = int(self.budget / max(1, 8 * itemsize * max(1, n_src)))
        return int(np.clip(rows, 16, 4096))

    def lam_at(self, tq_np, theta, W, src_hi=None):
        """lambda at query times, no grad. src_hi limits sources to self.t[:src_hi]."""
        mu, K, alpha, c, p = [float(x) for x in theta]
        n_src = self.n if src_hi is None else src_hi
        out = np.empty(len(tq_np))
        chunk = self._chunk(min(n_src, self._eff_sources(W)))
        with torch.no_grad():
            for x0 in range(0, len(tq_np), chunk):
                x1 = min(x0 + chunk, len(tq_np))
                tq = torch.as_tensor(tq_np[x0:x1], dtype=self.dtype, device=self.device)
                lo = int(np.searchsorted(self.t_np, tq_np[x0] - W, side="left")) if math.isfinite(W) else 0
                hi = n_src
                S = self._S_chunk_range(tq, alpha, c, p, W, lo, hi)
                out[x0:x1] = (mu + K * S).double().cpu().numpy()
        return out

    def _S_chunk_range(self, tq, alpha, c, p, W, lo, hi):
        a = alpha * LN10
        w = torch.exp(a * self.m[lo:hi])
        dt = tq[:, None] - self.t[None, lo:hi]
        valid = dt > 0
        if math.isfinite(W):
            valid = valid & (dt <= W)
        dtc = torch.where(valid, dt + c, torch.ones_like(dt))
        q = torch.exp(-p * torch.log(dtc)) * valid * w[None, :]
        return q.sum(dim=1)

    def _eff_sources(self, W):
        if not math.isfinite(W):
            return self.n
        span = self.t_np[-1] - self.t_np[0]
        return int(min(self.n, max(1024, self.n * W / max(span, 1e-9))))

    # -------- per-source integral of the kernel over [lo, hi], exp_h._G_terms in torch
    def _G(self, t_src, lo, hi, c, p, W):
        """t_src: (n_src,) torch; lo, hi: (rows,1) torch. Returns (rows, n_src)."""
        a_lo = torch.maximum(lo, t_src[None, :])
        a_hi = torch.minimum(hi, t_src[None, :] + W) if math.isfinite(W) else \
            hi.expand(-1, t_src.numel())
        live = a_hi > a_lo
        A = torch.where(live, a_hi - t_src[None, :] + c, torch.ones_like(a_lo))
        B = torch.where(live, a_lo - t_src[None, :] + c, torch.ones_like(a_lo))
        lA, lB = torch.log(A), torch.log(B)
        s = 1.0 - p
        if abs(float(s)) < 1e-6:
            G = (lA - lB) + s * (lA**2 - lB**2) / 2.0 + s * s * (lA**3 - lB**3) / 6.0
        else:
            G = (torch.exp(s * lA) - torch.exp(s * lB)) / s
        return G * live

    def seg_A(self, lo_np, hi_np, theta, W, src_hi=None):
        """Integral of lambda over each interval [lo_i, hi_i]. Returns numpy (n_int,)."""
        mu, K, alpha, c, p = [float(x) for x in theta]
        a = alpha * LN10
        n_src = self.n if src_hi is None else src_hi
        out = np.empty(len(lo_np))
        chunk = self._chunk(min(n_src, self._eff_sources(W)))
        with torch.no_grad():
            for x0 in range(0, len(lo_np), chunk):
                x1 = min(x0 + chunk, len(lo_np))
                lo = torch.as_tensor(lo_np[x0:x1], dtype=self.dtype, device=self.device)[:, None]
                hi = torch.as_tensor(hi_np[x0:x1], dtype=self.dtype, device=self.device)[:, None]
                # sources that can contribute: t_j < hi_max and (if truncated) t_j >= lo_min - W
                j0 = int(np.searchsorted(self.t_np, lo_np[x0] - W, side="left")) if math.isfinite(W) else 0
                j1 = min(n_src, int(np.searchsorted(self.t_np, hi_np[x1 - 1], side="right")))
                if j1 <= j0:
                    out[x0:x1] = mu * (hi_np[x0:x1] - lo_np[x0:x1])
                    continue
                w = torch.exp(a * self.m[j0:j1])
                G = self._G(self.t[j0:j1], lo, hi, c, p, W)
                out[x0:x1] = (mu * (hi - lo).squeeze(1) + K * (G * w[None, :]).sum(dim=1)).double().cpu().numpy()
        return out

    # -------- MLE on a target window
    def _neg_ll_and_grad(self, x, tgt_lo, tgt_hi, T0, T1, W):
        """-LL and d(-LL)/dx for x = log theta, targets self.t[tgt_lo:tgt_hi], window [T0,T1]."""
        xt = torch.as_tensor(x, dtype=self.dtype, device=self.device).requires_grad_(True)
        tq_np = self.t_np[tgt_lo:tgt_hi]
        chunk = self._chunk(self._eff_sources(W))
        total = 0.0
        bad = False
        for x0 in range(0, len(tq_np), chunk):
            x1 = min(x0 + chunk, len(tq_np))
            th = torch.exp(xt)                       # fresh graph per chunk: nothing retained
            mu, K, alpha, c, p = th[0], th[1], th[2], th[3], th[4]
            tq = self.t[tgt_lo + x0: tgt_lo + x1]
            lo = int(np.searchsorted(self.t_np, tq_np[x0] - W, side="left")) if math.isfinite(W) else 0
            S = self._S_chunk_range(tq, alpha, c, p, W, lo, tgt_lo + x1)
            lam = mu + K * S
            if not bool(torch.isfinite(lam).all()) or bool((lam <= 0).any()):
                bad = True
                break
            part = -torch.log(lam).sum()
            part.backward()
            total += float(part.detach().cpu())
            del S, lam, part, th
        if bad:
            return 1e18, np.zeros(5)
        # integral term: sources with t <= T1
        src_hi = int(np.searchsorted(self.t_np, T1, side="right"))
        lo_t = torch.as_tensor([T0], dtype=self.dtype, device=self.device)[:, None]
        hi_t = torch.as_tensor([T1], dtype=self.dtype, device=self.device)[:, None]
        step = max(1024, self._chunk(1))
        for j0 in range(0, src_hi, step):
            j1 = min(j0 + step, src_hi)
            th = torch.exp(xt)
            mu, K, alpha, c, p = th[0], th[1], th[2], th[3], th[4]
            w = torch.exp(alpha * LN10 * self.m[j0:j1])
            G = self._G(self.t[j0:j1], lo_t, hi_t, c, p, W)
            part = K * (G[0] * w).sum()
            if j0 == 0:
                part = part + mu * (T1 - T0)
            part.backward()
            total += float(part.detach().cpu())
            del G, w, part, th
        g = xt.grad.detach().double().cpu().numpy()
        v = float(total)
        if not np.isfinite(v) or not np.all(np.isfinite(g)):
            return 1e18, np.zeros(5)
        return v, g

    def fit(self, tgt_lo, tgt_hi, T0, T1, W, rate0, log=print, alpha_lb=None):
        lb = np.log([1e-6, 1e-6, ALPHA_BOUNDS[0] if alpha_lb is None else alpha_lb,
                     1e-6, P_BOUNDS[0]])
        ub = np.log([1e3, 1e2, ALPHA_BOUNDS[1], 1e1, P_BOUNDS[1]])
        bounds = list(zip(lb, ub))
        starts_out, best = [], None
        for si, (mu_f, K, al, c, p) in enumerate(STARTS):
            x0 = np.log(np.clip([mu_f * rate0, K, al, c, p], np.exp(lb) * 1.001, np.exp(ub) * 0.999))
            ts = time.time()
            r = minimize(lambda z: self._neg_ll_and_grad(z, tgt_lo, tgt_hi, T0, T1, W),
                         x0, jac=True, method="L-BFGS-B", bounds=bounds,
                         options={"maxiter": MAXITER, "maxfun": MAXITER * 2,
                                  "ftol": 1e-12, "gtol": 1e-8})
            th = np.exp(r.x)
            rec = {"start": dict(zip(["mu", "K", "alpha", "c", "p"],
                                     [float(mu_f * rate0), K, al, c, p])),
                   "LL": float(-r.fun),
                   "params": dict(zip(["mu", "K", "alpha", "c", "p"], map(float, th))),
                   "nit": int(r.nit), "success": bool(r.success), "seconds": round(time.time() - ts, 1)}
            starts_out.append(rec)
            log(f"[etas-fit] start {si+1}: LL={rec['LL']:.2f} mu={th[0]:.4f} K={th[1]:.5f} "
                f"alpha={th[2]:.3f} c={th[3]:.5f} p={th[4]:.4f} ({rec['nit']} it, {rec['seconds']}s)")
            if best is None or rec["LL"] > best["LL"]:
                best = rec
        theta = np.array([best["params"][k] for k in ["mu", "K", "alpha", "c", "p"]])
        return theta, starts_out, best


def verify_torch_matches_numpy(t, m, device, log=print):
    """Assert the torch ETAS value+gradient equals exp_h_etas.etas_ll on a subsample."""
    rng = np.random.default_rng(0)
    n = min(3000, len(t))
    ts = np.sort(t[:n])
    ms = m[:n]
    theta = np.array([0.5, 0.02, 1.0, 0.01, 1.15])
    T0, T1 = float(ts[500]), float(ts[-1])
    W = 200.0
    LL_np, g_np = eh.etas_ll(theta, ts, ms, 500, n, T0, T1, W, True, 512)
    E = TorchETAS(ts, ms, device)
    v, g_log = E._neg_ll_and_grad(np.log(theta), 500, n, T0, T1, W)
    LL_t = -v
    g_t = -g_log / theta          # d LL / d theta
    rel = abs(LL_t - LL_np) / max(1.0, abs(LL_np))
    grel = float(np.max(np.abs(g_t - g_np) / np.maximum(1.0, np.abs(g_np))))
    log(f"[verify] LL numpy={LL_np:.6f} torch={LL_t:.6f} rel={rel:.2e} | max grad rel={grel:.2e}")
    assert rel < 1e-9, f"torch/numpy ETAS LL mismatch: {rel}"
    assert grel < 1e-6, f"torch/numpy ETAS gradient mismatch: {grel}"
    return {"LL_numpy": float(LL_np), "LL_torch": float(LL_t), "rel_LL": float(rel),
            "max_rel_grad": grel}


# ----------------------------------------------------------------- solar hour + segments
def solar_hour(t_days, ref_lon):
    """Local solar hour-of-day (0..24) at a single reference longitude. t in days since epoch."""
    return np.mod((t_days + ref_lon / 360.0) * 24.0, 24.0)


def build_segments(t_bounds, ref_lon, n_bins=N_HOUR_BINS):
    """Split each interval [t_bounds[i], t_bounds[i+1]] at local-solar-hour-bin boundaries.

    Returns (lo, hi, bin_idx, interval_idx) as numpy arrays. sum over segments of an interval
    reconstructs the interval exactly.
    """
    width = 24.0 / n_bins / 24.0            # bin width in days
    off = ref_lon / 360.0                   # solar time = t + off (in days)
    los, his, bins, iidx = [], [], [], []
    for i in range(len(t_bounds) - 1):
        a, b = float(t_bounds[i]), float(t_bounds[i + 1])
        if b <= a:
            continue
        k0 = math.floor((a + off) / width)
        k1 = math.floor((b + off) / width)
        if k1 == k0:
            los.append(a); his.append(b); bins.append(k0 % n_bins); iidx.append(i)
        else:
            prev = a
            for k in range(k0, k1):
                edge = (k + 1) * width - off
                if edge > prev:
                    los.append(prev); his.append(edge); bins.append(k % n_bins); iidx.append(i)
                    prev = edge
            if b > prev:
                los.append(prev); his.append(b); bins.append(k1 % n_bins); iidx.append(i)
    return (np.array(los), np.array(his), np.array(bins, dtype=np.int64),
            np.array(iidx, dtype=np.int64))


# ----------------------------------------------------------------- neural TPP
class NeuralTPP(nn.Module):
    """GRU encoder -> log-normal mixture over tau_{i+1}. Strictly causal by construction."""

    def __init__(self, n_feat, hidden, layers, n_mix, dropout, mu_init=0.0, sd_init=1.0):
        super().__init__()
        self.n_mix = n_mix
        self.gru = nn.GRU(n_feat, hidden, num_layers=layers, batch_first=True,
                          dropout=dropout if layers > 1 else 0.0)
        self.head = nn.Linear(hidden, 3 * n_mix)
        nn.init.zeros_(self.head.weight)
        with torch.no_grad():
            b = torch.zeros(3 * n_mix)
            b[n_mix:2 * n_mix] = mu_init + sd_init * torch.linspace(-1.5, 1.5, n_mix)
            b[2 * n_mix:] = math.log(math.expm1(max(sd_init, 1e-2)))   # softplus^-1
            self.head.bias.copy_(b)

    def forward(self, x, h0=None):
        out, h = self.gru(x, h0)
        return self.head(out), h

    @staticmethod
    def log_f_tau(params, y, n_mix):
        """params (..., 3K); y = log(tau). Returns log density of TAU (per day)."""
        logits = params[..., :n_mix]
        mu = params[..., n_mix:2 * n_mix]
        sd = torch.nn.functional.softplus(params[..., 2 * n_mix:]) + 1e-4
        logw = torch.log_softmax(logits, dim=-1)
        z = (y[..., None] - mu) / sd
        logN = -0.5 * z * z - torch.log(sd) - 0.5 * math.log(2 * math.pi)
        return torch.logsumexp(logw + logN, dim=-1) - y      # jacobian d y / d tau = 1/tau


def make_features(t, mags, M0, ref_lon, use_mag=True, use_tod=True, norm=None):
    """Per-event input features. tau_i = t_i - t_{i-1} (tau_0 := median of the rest)."""
    tau = np.diff(t, prepend=t[0])
    pos = tau[tau > 0]
    tau[0] = np.median(pos) if len(pos) else 1.0
    tau = np.maximum(tau, TAU_FLOOR_DAYS)
    y = np.log(tau)
    mm = mags - M0
    h = solar_hour(t, ref_lon)
    if norm is None:
        norm = {"y_mean": float(y.mean()), "y_std": float(y.std() + 1e-12),
                "m_mean": float(mm.mean()), "m_std": float(mm.std() + 1e-12)}
    cols = [(y - norm["y_mean"]) / norm["y_std"]]
    if use_mag:
        cols.append((mm - norm["m_mean"]) / norm["m_std"])
    else:
        cols.append(np.zeros_like(mm))
    if use_tod:
        cols.append(np.sin(2 * np.pi * h / 24.0))
        cols.append(np.cos(2 * np.pi * h / 24.0))
    else:
        cols.append(np.zeros_like(h))
        cols.append(np.zeros_like(h))
    return np.stack(cols, axis=1).astype(np.float32), y.astype(np.float32), norm


def nn_full_pass(model, X, Y, device, n_mix, chunk=4096):
    """Left-to-right causal pass over the whole sequence, hidden state carried forward.

    Returns logf[i] = log f(tau_{i+1} | events <= i) placed at index i+1 (so logf[i] is the
    density of the observed tau_i given events < i). logf[0] is set to nan.
    """
    model.eval()
    N = len(X)
    out = np.full(N, np.nan)
    h = None
    Xt = torch.as_tensor(X, device=device)
    Yt = torch.as_tensor(Y, device=device)
    with torch.no_grad():
        for a in range(0, N - 1, chunk):
            b = min(a + chunk, N - 1)
            params, h = model(Xt[a:b][None, ...], h)
            lf = NeuralTPP.log_f_tau(params[0], Yt[a + 1:b + 1], n_mix)
            out[a + 1:b + 1] = lf.double().cpu().numpy()
    return out


def train_nn(X, Y, tr_lo, tr_hi, va_lo, va_hi, cfg, device, seed, log=print, tag=""):
    """Train with Adam, early stopping on validation LL. Returns (model, history)."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    n_mix = cfg["n_mix"]
    model = NeuralTPP(X.shape[1], cfg["hidden"], cfg["layers"], n_mix, cfg["dropout"],
                      mu_init=float(Y[tr_lo:tr_hi].mean()),
                      sd_init=float(Y[tr_lo:tr_hi].std() + 1e-3)).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=cfg["lr"])
    Xt = torch.as_tensor(X, device=device)
    Yt = torch.as_tensor(Y, device=device)
    L, Wm = cfg["chunk_len"], cfg["warmup_mask"]
    rng = np.random.default_rng(seed)
    hi_start = tr_hi - L                      # last legal chunk start: targets stay < tr_hi
    assert hi_start > 0, "training window too short for the chunk length"
    best = {"val": -np.inf, "step": -1, "state": None}
    hist, t0, stop = [], time.time(), None
    for step in range(1, cfg["max_steps"] + 1):
        model.train()
        starts = rng.integers(0, hi_start, size=cfg["batch"])
        idx = starts[:, None] + np.arange(L)[None, :]
        xb = Xt[torch.as_tensor(idx, device=device)]
        params, _ = model(xb)
        yb = Yt[torch.as_tensor(idx[:, 1:], device=device)]
        lf = NeuralTPP.log_f_tau(params[:, :-1, :], yb, n_mix)
        # mask: burn in Wm steps of context, and require the TARGET to lie in the train window
        tgt = idx[:, 1:]
        mask = np.zeros_like(tgt, dtype=bool)
        mask[:, Wm:] = True
        mask &= (tgt >= tr_lo) & (tgt < tr_hi)
        mt = torch.as_tensor(mask, device=device)
        if mt.sum() == 0:
            continue
        loss = -(lf * mt).sum() / mt.sum()
        opt.zero_grad(set_to_none=True)
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), cfg["grad_clip"])
        opt.step()
        if step % cfg["eval_every"] == 0 or step == cfg["max_steps"]:
            lf_all = nn_full_pass(model, X, Y, device, n_mix)
            val = float(np.mean(lf_all[va_lo:va_hi]))
            trn = float(np.mean(lf_all[tr_lo:tr_hi]))
            hist.append({"step": step, "train_logf": trn, "val_logf": val,
                         "minutes": round((time.time() - t0) / 60, 2)})
            log(f"[nn{tag}] step {step:>5d} train {trn:+.4f} val {val:+.4f} "
                f"({hist[-1]['minutes']:.1f} min)")
            if val > best["val"] + 1e-6:
                best = {"val": val, "step": step,
                        "state": {k: v.detach().clone() for k, v in model.state_dict().items()}}
            elif step - best["step"] >= cfg["patience"] * cfg["eval_every"]:
                stop = f"early stopping (patience {cfg['patience']} evals)"
                break
            if (time.time() - t0) / 60 > cfg["max_train_minutes"]:
                stop = f"time cap {cfg['max_train_minutes']} min"
                break
    if stop is None:
        stop = "max_steps reached"
    model.load_state_dict(best["state"])
    log(f"[nn{tag}] {stop}; best val {best['val']:+.4f} at step {best['step']}")
    return model, {"history": hist, "best_val_logf": best["val"], "best_step": best["step"],
                   "stop_reason": stop, "minutes": round((time.time() - t0) / 60, 2)}


# ----------------------------------------------------------------- ETAS simulator
def branching_ratio(theta, b, M0=None):
    mu, K, alpha, c, p = [float(x) for x in theta]
    if p <= 1 or b <= alpha:
        return None
    return K * c**(1 - p) / (p - 1) * (b / (b - alpha))


def simulate_poisson(rate, M0, b, T, n_target, t0_abs, rng):
    """Homogeneous Poisson + GR magnitudes. The DECISIVE leakage control: here the true
    inter-event density is exponential and independent of history, so a causal model can earn
    no bits over a correctly-fitted ETAS (which should collapse to K ~ 0)."""
    n = n_target
    t = np.sort(rng.uniform(0, T, size=n))
    m = M0 + rng.exponential(1.0 / (b * LN10), size=n)
    return t + t0_abs, m, {"generator": "homogeneous Poisson", "rate_per_day": float(n / T),
                           "n_total": int(n)}


def simulate_etas(theta, M0, b, T, n_target, t0_abs, rng, max_events=None):
    """Branching-process ETAS simulation. mu is rescaled so E[n] ~= n_target."""
    mu, K, alpha, c, p = [float(x) for x in theta]
    n_br = K * c**(1 - p) / (p - 1) * (b / (b - alpha)) if (p > 1 and b > alpha) else None
    if n_br is not None and 0 < n_br < 0.98:
        mu_sim = n_target * (1 - n_br) / T
    else:
        mu_sim = mu
        n_br = None
    max_events = max_events or int(4 * n_target)

    def gr_mags(n):
        return M0 + rng.exponential(1.0 / (b * LN10), size=n)

    n0 = rng.poisson(mu_sim * T)
    times = list(np.sort(rng.uniform(0, T, size=n0)))
    mags = list(gr_mags(n0))
    cur_t = np.array(times)
    cur_m = np.array(mags)
    gen = 0
    while len(cur_t) and len(times) < max_events:
        gen += 1
        w = K * 10 ** (alpha * (cur_m - M0))
        rem = T - cur_t
        expn = w * ((rem + c) ** (1 - p) - c ** (1 - p)) / (1 - p)
        nk = rng.poisson(np.maximum(expn, 0))
        tot = int(nk.sum())
        if tot == 0:
            break
        parent = np.repeat(np.arange(len(cur_t)), nk)
        u = rng.uniform(size=tot)
        Cp = c ** (1 - p)
        Rp = (rem[parent] + c) ** (1 - p)
        dt = (Cp + u * (Rp - Cp)) ** (1 / (1 - p)) - c
        ct = cur_t[parent] + dt
        keep = (ct > 0) & (ct < T)
        ct = ct[keep]
        cm = gr_mags(len(ct))
        times.extend(ct.tolist())
        mags.extend(cm.tolist())
        cur_t, cur_m = ct, cm
        if gen > 200:
            break
    t = np.array(times)
    m = np.array(mags)
    o = np.argsort(t)
    t, m = t[o], m[o]
    info = {"mu_sim": float(mu_sim), "branching_ratio_used": n_br, "generations": gen,
            "n_background": int(n0), "n_simulated": int(len(t))}
    if len(t) > n_target:
        # a supercritical realisation overshoots; cut the TIME WINDOW at the n_target-th
        # event so the control catalogue has the intended size. Still exactly an ETAS
        # realisation, on a shorter window.
        cut = float(t[n_target - 1])
        t, m = t[:n_target], m[:n_target]
        info.update(time_truncated_to_days=cut, time_truncated=True)
    else:
        info["time_truncated"] = False
    info["n_total"] = int(len(t))
    return t + t0_abs, m, info


# ----------------------------------------------------------------- one full pipeline
def run_pipeline(name, t_all, m_all, M0, ref_lon, device, dev_info, cfg, seed,
                 do_ablations=False, log=print, alpha_lb=None):
    """Fit + score every model on one catalogue. Returns a results dict."""
    R = {"name": name, "n_events": int(len(t_all)), "M0": M0,
         "reference_longitude_deg": ref_lon}
    t0_run = time.time()

    # ---- splits (by TIME fraction, engine/splits.py rule) ----
    T_first, T_last = float(t_all[0]), float(t_all[-1])
    span = T_last - T_first
    T_expl = T_first + EXPLORE_FRAC * span
    e_span = T_expl - T_first
    T_tr = T_first + TRAIN_FRAC * e_span
    T_va = T_first + (TRAIN_FRAC + VAL_FRAC) * e_span
    tr_hi = int(np.searchsorted(t_all, T_tr, side="left"))
    va_hi = int(np.searchsorted(t_all, T_va, side="left"))
    te_hi = int(np.searchsorted(t_all, T_expl, side="left"))
    T0_fit = T_first + BURN_IN_DAYS
    fit_lo = int(np.searchsorted(t_all, T0_fit, side="left"))
    R["split"] = {
        "rule": "exploration = first 70% of span; within it train 60% / val 10% / test 30% "
                "by TIME fraction of the exploration span",
        "t_first": T_first, "t_explore_end": T_expl, "t_train_end": T_tr, "t_val_end": T_va,
        "n_train": tr_hi, "n_val": va_hi - tr_hi, "n_test": te_hi - va_hi,
        "n_exploration": te_hi, "n_holdout_never_read": int(len(t_all) - te_hi),
        "etas_burn_in_days": BURN_IN_DAYS, "n_etas_fit_targets": tr_hi - fit_lo,
    }
    # ---- counted invariants ----
    assert 0 < fit_lo < tr_hi < va_hi < te_hi <= len(t_all), "degenerate split"
    assert tr_hi + (va_hi - tr_hi) + (te_hi - va_hi) == te_hi, "split counts do not sum"
    assert t_all[va_hi] > t_all[tr_hi - 1], "TEST timestamp <= max TRAIN timestamp"
    assert float(t_all[va_hi]) > float(t_all[:tr_hi].max()), "future leakage in the split"
    assert np.all(np.diff(t_all) >= 0), "catalogue not sorted in time"
    log(f"[{name}] n={len(t_all)} | train {tr_hi} (fit targets {tr_hi-fit_lo}) | "
        f"val {va_hi-tr_hi} | test {te_hi-va_hi} | holdout {len(t_all)-te_hi} (never read)")

    # ---- ETAS fit on train ----
    E = TorchETAS(t_all[:te_hi], m_all[:te_hi] - M0, device)
    rate_fit = (tr_hi - fit_lo) / (T_tr - T0_fit)
    ts = time.time()
    theta, starts_out, best = E.fit(fit_lo, tr_hi, T0_fit, T_tr, W_FIT_DAYS, rate_fit,
                                    log=log, alpha_lb=alpha_lb)
    mu, K, alpha, c, p = [float(x) for x in theta]
    b_train = eh.aki_b(m_all[:tr_hi], M0)
    n_br = (K * c**(1 - p) / (p - 1) * (b_train / (b_train - alpha))
            if (p > 1 and b_train > alpha) else None)
    R["etas_fit"] = {"alpha_lower_bound": ALPHA_BOUNDS[0] if alpha_lb is None else alpha_lb,
                     "alpha_at_lower_bound": bool(abs(alpha - (ALPHA_BOUNDS[0] if alpha_lb is None else alpha_lb)) < 1e-6),
                     "frozen_params": {"mu": mu, "K": K, "alpha": alpha, "c": c, "p": p, "M0": M0},
                     "truncation_W_days_fit": W_FIT_DAYS, "train_rate_per_day": float(rate_fit),
                     "b_value_train_aki": b_train, "branching_ratio_n": n_br,
                     "starts": starts_out, "best_LL": best["LL"],
                     "seconds": round(time.time() - ts, 1)}
    log(f"[{name}] ETAS FROZEN mu={mu:.4f} K={K:.5f} alpha={alpha:.3f} c={c:.5f} p={p:.4f} "
        f"| b={b_train:.3f} n={n_br}")

    # ---- segments (exact per-interval, hour-bin-split ETAS integrals) ----
    def segment_scores(lo_idx, hi_idx, W):
        """Per-event ETAS log f(tau_i) and per-bin integrals for events lo_idx..hi_idx-1."""
        bounds = t_all[lo_idx - 1: hi_idx]                    # tiling of the window
        slo, shi, sbin, sidx = build_segments(bounds, ref_lon)
        A = E.seg_A(slo, shi, theta, W, src_hi=hi_idx)
        assert np.all(np.isfinite(A)), "non-finite ETAS interval integral"
        n_int = hi_idx - lo_idx
        Lam = np.zeros(n_int)
        np.add.at(Lam, sidx, A)
        lam = E.lam_at(t_all[lo_idx:hi_idx], theta, W, src_hi=hi_idx)
        assert np.all(lam > 0) and np.all(np.isfinite(lam)), "non-positive/non-finite lambda"
        return lam, Lam, (slo, shi, sbin, sidx, A)

    # train: fit the diurnal factor by its exact MLE g_b = n_b / I_b
    ts = time.time()
    lam_tr, Lam_tr, seg_tr = segment_scores(fit_lo, tr_hi, W_SCORE_DAYS)
    slo, shi, sbin, sidx, A_tr = seg_tr
    I_b = np.zeros(N_HOUR_BINS)
    np.add.at(I_b, sbin, A_tr)
    h_ev_tr = solar_hour(t_all[fit_lo:tr_hi], ref_lon)
    bin_ev_tr = np.floor(h_ev_tr / (24.0 / N_HOUR_BINS)).astype(int) % N_HOUR_BINS
    n_b = np.bincount(bin_ev_tr, minlength=N_HOUR_BINS).astype(float)
    assert np.all(I_b > 0), "empty hour-bin integral"
    g = n_b / I_b
    R["diurnal"] = {
        "n_bins": N_HOUR_BINS,
        "fit": "exact MLE g_b = n_b / I_b on the train window, same truncation W as scoring",
        "counts": n_b.tolist(), "etas_integral_per_bin": I_b.tolist(), "g": g.tolist(),
        "g_normalised": (g / g.mean()).tolist(),
        "amplitude_max_over_min": float(g.max() / g.min()),
        "train_LL_gain_bits_per_event": float(
            (np.sum(np.log(g[bin_ev_tr])) - np.sum((g - 1.0) * I_b)) / (tr_hi - fit_lo) / LN2),
        "seconds": round(time.time() - ts, 1),
    }
    log(f"[{name}] diurnal g: max/min={g.max()/g.min():.3f}, "
        f"train gain {R['diurnal']['train_LL_gain_bits_per_event']:+.4f} bits/event")

    # ---- test scoring, walk-forward, history = all prior events ----
    ts = time.time()
    lam_te, Lam_te, seg_te = segment_scores(va_hi, te_hi, W_SCORE_DAYS)
    slo, shi, sbin, sidx, A_te = seg_te
    n_te = te_hi - va_hi
    tau_te = np.diff(t_all[va_hi - 1: te_hi])
    assert np.all(tau_te >= 0), "negative inter-event time"

    logf_etas = np.log(lam_te) - Lam_te
    # diurnal: log lambda gets log g at the event's bin; the integral uses g per segment
    h_ev_te = solar_hour(t_all[va_hi:te_hi], ref_lon)
    bin_ev_te = np.floor(h_ev_te / (24.0 / N_HOUR_BINS)).astype(int) % N_HOUR_BINS
    Lam_te_d = np.zeros(n_te)
    np.add.at(Lam_te_d, sidx, A_te * g[sbin])
    logf_diurnal = np.log(lam_te * g[bin_ev_te]) - Lam_te_d
    rate_poisson = float(rate_fit)
    logf_poisson = np.log(rate_poisson) - rate_poisson * tau_te
    for nm, v in [("etas", logf_etas), ("diurnal", logf_diurnal), ("poisson", logf_poisson)]:
        assert np.all(np.isfinite(v)), f"non-finite log f for {nm}"
    R["test_scoring"] = {
        "n_test_events": n_te, "window_days": float(t_all[te_hi - 1] - t_all[va_hi - 1]),
        "truncation_W_days_score": W_SCORE_DAYS,
        "definition": "log f(tau_i) = log lambda(t_i) - Integral_{t_{i-1}}^{t_i} lambda; "
                      "intervals tile the window exactly, so bits/event is a plain mean",
        "seconds": round(time.time() - ts, 1),
    }
    # truncation cost: re-score a random subsample untruncated
    rng = np.random.default_rng(seed)
    sub = np.sort(rng.choice(n_te, size=min(2000, n_te), replace=False))
    lam_u = E.lam_at(t_all[va_hi:te_hi][sub], theta, np.inf, src_hi=te_hi)
    lo_u = t_all[va_hi - 1: te_hi - 1][sub]
    hi_u = t_all[va_hi:te_hi][sub]
    Lam_u = E.seg_A(lo_u, hi_u, theta, np.inf, src_hi=te_hi)
    d_trunc = float(np.mean((np.log(lam_u) - Lam_u) - logf_etas[sub]) / LN2)
    R["test_scoring"]["truncation_cost_bits_per_event"] = d_trunc
    R["test_scoring"]["truncation_note"] = (
        "untruncated minus truncated ETAS log f on a 2000-event subsample; positive means the "
        "truncated ETAS baseline is weakened by this many bits/event")
    log(f"[{name}] truncation cost {d_trunc:+.5f} bits/event (W={W_SCORE_DAYS} d)")

    # ---- neural TPP ----
    def nn_run(use_mag, use_tod, tag, sd):
        # standardisation constants must come from TRAIN ONLY
        _, _, norm_tr = make_features(t_all[:tr_hi], m_all[:tr_hi], M0, ref_lon,
                                      use_mag=use_mag, use_tod=use_tod, norm=None)
        X, Y, _ = make_features(t_all[:te_hi], m_all[:te_hi], M0, ref_lon,
                                use_mag=use_mag, use_tod=use_tod, norm=norm_tr)
        model, hist = train_nn(X, Y, fit_lo, tr_hi, tr_hi, va_hi, cfg, device, sd,
                               log=log, tag=tag)
        lf = nn_full_pass(model, X, Y, device, cfg["n_mix"])
        assert np.all(np.isfinite(lf[1:])), "non-finite neural log f"
        return model, hist, lf[va_hi:te_hi], norm_tr

    model, nn_hist, logf_nn, norm_tr = nn_run(True, True, "", seed)

    # tau == 0 (tied timestamps) has infinite density under any continuous model: drop for ALL
    ok = tau_te > 0
    n_dropped = int((~ok).sum())
    logf_etas, logf_diurnal, logf_poisson, logf_nn = (
        logf_etas[ok], logf_diurnal[ok], logf_poisson[ok], logf_nn[ok])
    t_te = t_all[va_hi:te_hi][ok]
    n_te_scored = int(ok.sum())
    assert n_te_scored == n_te - n_dropped
    R["test_scoring"].update(
        n_tied_timestamps_dropped=n_dropped, n_scored=n_te_scored,
        tie_rule="events with tau == 0 are dropped from scoring for EVERY model",
        mean_logf={"etas": float(logf_etas.mean()), "etas_diurnal": float(logf_diurnal.mean()),
                   "poisson": float(logf_poisson.mean()), "neural": float(logf_nn.mean())})
    R["neural"] = {"config": {k: v for k, v in cfg.items()},
                   "n_features": 4, "inputs": ["z(log tau)", "z(M-M0)", "sin(2pi h/24)", "cos(2pi h/24)"],
                   "normalisation_train_only": norm_tr,
                   "training": nn_hist, "seed": seed,
                   "n_parameters": int(sum(q.numel() for q in model.parameters()))}

    # ---- the numbers ----
    def bits(a, b_):
        return float(np.mean(a - b_) / LN2)

    mags_te = m_all[va_hi:te_hi][ok]
    big = mags_te >= 2.5 - 1e-9
    R["bits_per_event"] = {
        "neural_over_etas_diurnal_PRIMARY": bits(logf_nn, logf_diurnal),
        "neural_over_etas": bits(logf_nn, logf_etas),
        "neural_over_poisson": bits(logf_nn, logf_poisson),
        "etas_diurnal_over_etas": bits(logf_diurnal, logf_etas),
        "etas_over_poisson": bits(logf_etas, logf_poisson),
        "n_test_events": n_te_scored,
        "M>=2.5_subset": {
            "n": int(big.sum()),
            "neural_over_etas_diurnal": bits(logf_nn[big], logf_diurnal[big]) if big.sum() else None,
            "neural_over_etas": bits(logf_nn[big], logf_etas[big]) if big.sum() else None,
            "etas_diurnal_over_etas": bits(logf_diurnal[big], logf_etas[big]) if big.sum() else None,
        },
    }
    # a crude uncertainty: block bootstrap over 200 contiguous blocks of test events
    d = (logf_nn - logf_diurnal) / LN2
    nb = 200
    edges = np.linspace(0, n_te_scored, nb + 1).astype(int)
    blocks = [d[edges[i]:edges[i + 1]].sum() for i in range(nb)]
    sizes = [edges[i + 1] - edges[i] for i in range(nb)]
    boot = []
    for _ in range(2000):
        pick = rng.integers(0, nb, size=nb)
        boot.append(np.sum([blocks[j] for j in pick]) / np.sum([sizes[j] for j in pick]))
    R["bits_per_event"]["primary_block_bootstrap_95CI"] = [float(np.percentile(boot, 2.5)),
                                                           float(np.percentile(boot, 97.5))]
    R["bits_per_event"]["bootstrap_note"] = ("200 contiguous blocks of test events, 2000 resamples; "
                                             "descriptive, not a gate")

    log(f"[{name}] PRIMARY neural over ETAS+diurnal = "
        f"{R['bits_per_event']['neural_over_etas_diurnal_PRIMARY']:+.4f} bits/event")
    log(f"[{name}] neural over raw ETAS            = {R['bits_per_event']['neural_over_etas']:+.4f}")
    log(f"[{name}] ETAS+diurnal over ETAS          = {R['bits_per_event']['etas_diurnal_over_etas']:+.4f}")
    log(f"[{name}] neural over Poisson             = {R['bits_per_event']['neural_over_poisson']:+.4f}")
    log(f"[{name}] M>=2.5 (n={int(big.sum())}) neural over ETAS+diurnal = "
        f"{R['bits_per_event']['M>=2.5_subset']['neural_over_etas_diurnal']}")

    # ---- exploratory attribution (only if the neural model wins) ----
    won = R["bits_per_event"]["neural_over_etas_diurnal_PRIMARY"] > SUCCESS_BITS
    R["attribution_run"] = bool(won and do_ablations)
    if won and do_ablations:
        attr = {"NOTE": "EXPLORATORY. Not pre-registered, not multiplicity-corrected, "
                        "claims nothing."}
        bands = [("<1.5", 0.0, 1.5), ("1.5-2", 1.5, 2.0), ("2-2.5", 2.0, 2.5),
                 ("2.5-3", 2.5, 3.0), (">=3", 3.0, 99.0)]
        attr["by_magnitude_band"] = {}
        for bn, a_, b_ in bands:
            s = (mags_te >= a_ - 1e-9) & (mags_te < b_)
            attr["by_magnitude_band"][bn] = {
                "n": int(s.sum()),
                "over_etas_diurnal": bits(logf_nn[s], logf_diurnal[s]) if s.sum() > 20 else None,
                "over_etas": bits(logf_nn[s], logf_etas[s]) if s.sum() > 20 else None}
        # time since the last M>=4 anywhere in the catalogue (history includes train)
        big4 = np.where(m_all[:te_hi] >= 4.0 - 1e-9)[0]
        if len(big4):
            j = np.searchsorted(t_all[big4], t_te, side="right") - 1
            has = j >= 0
            dt4 = np.full(len(t_te), np.inf)
            dt4[has] = t_te[has] - t_all[big4][j[has]]
            edges4 = [0, 1, 7, 30, 180, 1e9]
            attr["by_days_since_last_M4"] = {}
            for i in range(len(edges4) - 1):
                s = (dt4 >= edges4[i]) & (dt4 < edges4[i + 1])
                attr["by_days_since_last_M4"][f"{edges4[i]}-{edges4[i+1]}d"] = {
                    "n": int(s.sum()),
                    "over_etas_diurnal": bits(logf_nn[s], logf_diurnal[s]) if s.sum() > 20 else None}
        # ablations
        for tag, um, ut in [("_no_mag", False, True), ("_no_tod", True, False)]:
            _, h2, lf2, _ = nn_run(um, ut, tag, seed + 1)
            lf2 = lf2[ok]
            attr["ablation" + tag] = {
                "over_etas_diurnal": bits(lf2, logf_diurnal),
                "over_etas": bits(lf2, logf_etas),
                "over_full_model": bits(lf2, logf_nn),
                "best_val_logf": h2["best_val_logf"], "stop_reason": h2["stop_reason"]}
            log(f"[{name}] ablation{tag}: over ETAS+diurnal "
                f"{attr['ablation'+tag]['over_etas_diurnal']:+.4f}")
        R["attribution_exploratory"] = attr

    R["runtime_minutes"] = round((time.time() - t0_run) / 60, 2)
    R["_model"] = model
    R["_arrays"] = dict(mags_te=mags_te, t_te=t_te, logf_nn=logf_nn, logf_etas=logf_etas, logf_diurnal=logf_diurnal,
                        logf_poisson=logf_poisson, theta=theta, b_train=b_train,
                        T_span_train=(T0_fit, T_tr), n_all=len(t_all))
    return R


# ----------------------------------------------------------------- main
def main():
    t_start = time.time()
    MODELDIR.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    device, dev_info = pick_device()
    cfg = dict(NN) if device.type == "cuda" else dict(NN_CPU)
    print(f"[device] {dev_info['chosen']} :: {dev_info['reason']}")
    for d in dev_info["devices"]:
        print(f"[device] cuda:{d['index']} {d['name']} free {d['free_MiB']} MiB")

    res = {"experiment": "ARM-C neural TPP vs ETAS",
           "state_class": "first-run, exploration split",
           "run_utc": pd.Timestamp.utcnow().isoformat(),
           "seed": SEED, "torch_version": torch.__version__,
           "device": dev_info,
           "gpu_used": device.type == "cuda",
           "frozen": {
               "split": "exploration = first 70% of span; train 60% / val 10% / test 30% of "
                        "the exploration span, by time",
               "primary_statistic": "bits/event of the neural TPP over ETAS+diurnal on the "
                                    "QTM test window, time-only per-event log-densities",
               "success_rule": f"> {SUCCESS_BITS} bits/event on QTM",
               "calibration_control_rule": f"<= ~{CONTROL_TOL_BITS} bits/event on ETAS-simulated data",
               "etas_form": "exp_h_etas.py: lambda = mu + sum K 10^{alpha (M-M0)} (t-ti+c)^{-p}",
               "nn_config_gpu": NN, "nn_config_cpu_fallback": NN_CPU,
           },
           "catalogs": {}, "notes": []}

    # ------------------------------------------------ QTM (primary)
    spec = CATALOGS["QTM"]
    df, order = eh.load_catalog(spec["file"])
    d = df[df.mag >= spec["M0"] - 1e-9].sort_values("t").reset_index(drop=True)
    t_all = d.t.to_numpy(float)
    m_all = d.mag.to_numpy(float)
    ref_lon = float(d.lon.median())
    res["catalogs"]["QTM"] = {"file": spec["file"], "detected_column_order": order,
                              "n_rows_raw": int(len(df)), "n_at_floor": int(len(d)),
                              "M0": spec["M0"], "t_first": str(d.ts.iloc[0]),
                              "t_last": str(d.ts.iloc[-1]), "reference_longitude": ref_lon}
    print(f"[QTM] {len(df)} rows ({order} order), {len(d)} at M>={spec['M0']}, ref lon {ref_lon:.3f}")

    ver = verify_torch_matches_numpy(t_all[:6000], m_all[:6000] - spec["M0"], device)
    res["torch_vs_numpy_etas_check"] = ver

    # ------------------------------------------------ calibration control FIRST
    print("\n" + "=" * 70 + "\n[control] ETAS-simulated catalogue, identical pipeline\n" + "=" * 70)
    # a quick ETAS fit on the real train window is needed to seed the simulator; do the real
    # run first so the fitted parameters exist, then simulate from them.
    try:
        qtm = run_pipeline("QTM", t_all, m_all, spec["M0"], ref_lon, device, dev_info,
                           cfg, SEED, do_ablations=True)
    except torch.cuda.OutOfMemoryError as exc:
        res["notes"].append(f"GPU OOM on QTM ({exc}); fell back to CPU with the reduced config")
        torch.cuda.empty_cache()
        device = torch.device("cpu")
        cfg = dict(NN_CPU)
        res["device"]["chosen"] = "cpu (OOM fallback)"
        res["gpu_used"] = False
        qtm = run_pipeline("QTM", t_all, m_all, spec["M0"], ref_lon, device, dev_info,
                           cfg, SEED, do_ablations=True)
    torch.save(qtm.pop("_model").state_dict(), MODELDIR / "qtm_neural_tpp.pt")
    arrays = qtm.pop("_arrays")
    np.savez_compressed(MODELDIR / "qtm_test_logf.npz",
                        **{k: v for k, v in arrays.items() if isinstance(v, np.ndarray)})
    res["catalogs"]["QTM"].update(qtm)

    # ---- calibration controls: three generators ----
    theta = arrays["theta"]
    b_tr = arrays["b_train"]
    T_sim = float(t_all[-1] - t_all[0])
    n_br_real = branching_ratio(theta, b_tr)
    ctrl = {"rule": f"on data generated by a KNOWN model the neural TPP must not beat ETAS by "
                    f"more than ~{CONTROL_TOL_BITS} bits/event",
            "fitted_branching_ratio_on_real_QTM": n_br_real,
            "note": ("the QTM M>=1.0 ETAS fit is SUPERCRITICAL (n>1), so a same-size stationary "
                     "ETAS simulation from the fitted parameters is impossible; three "
                     "generators are used instead"),
            "generators": {}}

    def run_control(label, ts_, ms_, siminfo, sd):
        rec = {"sim": siminfo}
        try:
            cres = run_pipeline(label, ts_, ms_, spec["M0"], ref_lon, device, dev_info,
                                cfg, sd, do_ablations=False)
            cres.pop("_model"); cres.pop("_arrays")
            bp = cres["bits_per_event"]
            rec["neural_over_etas_bits_per_event"] = bp["neural_over_etas"]
            rec["neural_over_etas_diurnal_bits_per_event"] = bp["neural_over_etas_diurnal_PRIMARY"]
            rec["etas_over_poisson_bits_per_event"] = bp["etas_over_poisson"]
            rec["etas_fit"] = cres["etas_fit"]["frozen_params"]
            rec["n_test_scored"] = cres["test_scoring"]["n_scored"]
            rec["PASS"] = bool(bp["neural_over_etas"] <= CONTROL_TOL_BITS)
            print(f"[control:{label}] neural over ETAS = {bp['neural_over_etas']:+.4f} "
                  f"bits/event -> {'PASS' if rec['PASS'] else 'FAIL'}")
        except Exception as exc:  # noqa: BLE001
            rec["error"] = f"{type(exc).__name__}: {exc}"
            rec["PASS"] = None
            print(f"[control:{label}] FAILED TO RUN: {rec['error']}")
        return rec

    # (C) homogeneous Poisson -- the decisive leakage test
    rngc = np.random.default_rng(SEED + 11)
    ts_, ms_, si = simulate_poisson(len(t_all) / T_sim, spec["M0"], b_tr, T_sim, len(t_all),
                                    float(t_all[0]), rngc)
    ctrl["generators"]["poisson"] = run_control("CTRL-poisson", ts_, ms_, si, SEED + 11)

    # (B) SUBCRITICAL ETAS: K rescaled to branching ratio 0.90, same window and size
    if n_br_real and n_br_real > 0:
        th_b = np.array(theta, float).copy()
        th_b[1] = theta[1] * 0.90 / n_br_real
        rngc = np.random.default_rng(SEED + 13)
        ts_, ms_, si = simulate_etas(th_b, spec["M0"], b_tr, T_sim, len(t_all),
                                     float(t_all[0]), rngc)
        si["params"] = dict(zip(["mu", "K", "alpha", "c", "p"], map(float, th_b)))
        si["rescale"] = "K scaled so the branching ratio is 0.90 (the fitted value is >1)"
        ctrl["generators"]["etas_subcritical_n0.90"] = run_control(
            "CTRL-etas-sub", ts_, ms_, si, SEED + 13)

    # (A) ETAS at the fitted parameters, time-truncated to the target size
    rngc = np.random.default_rng(SEED + 7)
    ts_, ms_, si = simulate_etas(theta, spec["M0"], b_tr, T_sim, len(t_all),
                                 float(t_all[0]), rngc)
    si["params"] = dict(zip(["mu", "K", "alpha", "c", "p"], map(float, theta)))
    ctrl["generators"]["etas_fitted_supercritical_timecut"] = run_control(
        "CTRL-etas-fit", ts_, ms_, si, SEED + 7)

    passes = [v.get("PASS") for v in ctrl["generators"].values()]
    ctrl["PASS"] = bool(all(x is True for x in passes))
    ctrl["PASS_poisson_leakage_test"] = ctrl["generators"]["poisson"].get("PASS")
    ctrl["neural_over_etas_bits_per_event"] = {
        k: v.get("neural_over_etas_bits_per_event") for k, v in ctrl["generators"].items()}
    res["calibration_control"] = ctrl

    # ------------------------------------------------ ETAS alpha-bound sensitivity
    # The exp_h bounds pin alpha at its lower limit of 0.5 on QTM at M>=1.0. Relaxing the
    # bound gives ETAS a STRONGER fit, which can only shrink the neural gain -- the
    # conservative direction. Reported alongside, never instead of, the frozen number.
    try:
        alt = run_pipeline("QTM-alpha-relaxed", t_all, m_all, spec["M0"], ref_lon, device,
                           dev_info, cfg, SEED, do_ablations=False, alpha_lb=0.05)
        alt.pop("_model"); alt.pop("_arrays")
        res["etas_alpha_sensitivity"] = {
            "why": "the exp_h alpha lower bound (0.5) is ACTIVE on QTM at M>=1.0; refitting "
                   "with a lower bound of 0.05 gives ETAS a strictly better train fit",
            "alpha_lower_bound": 0.05,
            "etas_params": alt["etas_fit"]["frozen_params"],
            "alpha_at_lower_bound": alt["etas_fit"]["alpha_at_lower_bound"],
            "neural_over_etas_diurnal": alt["bits_per_event"]["neural_over_etas_diurnal_PRIMARY"],
            "neural_over_etas": alt["bits_per_event"]["neural_over_etas"],
            "etas_over_poisson": alt["bits_per_event"]["etas_over_poisson"],
            "M>=2.5": alt["bits_per_event"]["M>=2.5_subset"]}
        print(f"[alpha-sens] neural over ETAS+diurnal with alpha_lb=0.05: "
              f"{alt['bits_per_event']['neural_over_etas_diurnal_PRIMARY']:+.4f} bits/event")
    except Exception as exc:  # noqa: BLE001
        res["etas_alpha_sensitivity"] = {"error": f"{type(exc).__name__}: {exc}"}
        print(f"[alpha-sens] failed: {exc}")

    # ------------------------------------------------ verdict
    prim = res["catalogs"]["QTM"]["bits_per_event"]["neural_over_etas_diurnal_PRIMARY"]
    res["verdict"] = {
        "primary_bits_per_event_over_etas_diurnal": prim,
        "over_raw_etas": res["catalogs"]["QTM"]["bits_per_event"]["neural_over_etas"],
        "success_rule": f"> {SUCCESS_BITS} bits/event over ETAS+diurnal on QTM",
        "PASS": bool(prim > SUCCESS_BITS),
        "calibration_control_PASS": ctrl.get("PASS"),
        "reportable": bool(prim > SUCCESS_BITS and ctrl.get("PASS") is True),
    }
    res["runtime_minutes_total"] = round((time.time() - t_start) / 60, 2)
    OUT.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    print("\n" + "=" * 70)
    print(f"PRIMARY (QTM, neural over ETAS+diurnal): {prim:+.4f} bits/event  "
          f"[{'PASS' if prim > SUCCESS_BITS else 'FAIL'} vs {SUCCESS_BITS}]")
    for k, v in ctrl.get("neural_over_etas_bits_per_event", {}).items():
        print(f"  control [{k}]: neural over ETAS {v:+.4f} bits/event")
    print(f"calibration control overall: {ctrl.get('PASS')} "
          f"(Poisson leakage test: {ctrl.get('PASS_poisson_leakage_test')})")
    print(f"runtime {res['runtime_minutes_total']} min -> {OUT.name}")
    print("=" * 70)

    # ------------------------------------------------ SCSN (secondary, if time)
    if (time.time() - t_start) / 60 < 170:
        try:
            spec = CATALOGS["SCSN"]
            df, order = eh.load_catalog(spec["file"])
            box = df.lat.between(*spec["box"]["lat"]) & df.lon.between(*spec["box"]["lon"])
            d = df[box & (df.mag >= spec["M0"] - 1e-9)].sort_values("t").reset_index(drop=True)
            ref_lon2 = float(d.lon.median())
            print(f"\n[SCSN] {len(d)} events at M>={spec['M0']} in the exp_h box")
            # attribution is enabled here because SCSN is real data; it runs only if the
            # neural model actually clears the gate, and it is EXPLORATORY either way
            s = run_pipeline("SCSN", d.t.to_numpy(float), d.mag.to_numpy(float), spec["M0"],
                             ref_lon2, device, dev_info, cfg, SEED, do_ablations=True)
            torch.save(s.pop("_model").state_dict(), MODELDIR / "scsn_neural_tpp.pt")
            s.pop("_arrays")
            res["catalogs"]["SCSN"] = {"file": spec["file"], "n_at_floor": int(len(d)),
                                       "box": spec["box"], "reference_longitude": ref_lon2, **s}
        except Exception as exc:  # noqa: BLE001
            res["catalogs"]["SCSN"] = {"error": f"{type(exc).__name__}: {exc}"}
            print(f"[SCSN] failed: {exc}")
    else:
        res["catalogs"]["SCSN"] = {"skipped": "time budget"}

    res["runtime_minutes_total"] = round((time.time() - t_start) / 60, 2)
    OUT.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    print(f"\n[done] total runtime {res['runtime_minutes_total']} min -> {OUT.name}")


if __name__ == "__main__":
    main()
