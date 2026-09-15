import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch


def nb_nll(x, mu, phi):
    size = 1.0 / phi
    return -(
        torch.lgamma(x + size)
        - torch.lgamma(size)
        - torch.lgamma(x + 1.0)
        + size * (torch.log(size) - torch.log(size + mu))
        + x * (torch.log(mu) - torch.log(size + mu))
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--mu0", required=True)
    ap.add_argument("--phi", required=True)
    ap.add_argument("--clip", type=float, required=True)
    ap.add_argument("--conditions", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rank", type=int, default=10)
    ap.add_argument("--steps", type=int, default=1000)
    ap.add_argument("--batch-size", type=int, default=4096)
    ap.add_argument("--seed", type=int, default=20260829)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    x_np = np.load(args.counts).astype(np.float32)
    mu0_np = np.load(args.mu0).astype(np.float32)
    phi_np = np.load(args.phi).astype(np.float32).reshape(-1)
    meta = pd.read_csv(args.meta, dtype=str)
    conditions = args.conditions.split(",")
    cond = pd.Categorical(meta["condition"], categories=conditions).codes
    if np.any(cond < 0) or x_np.shape != mu0_np.shape:
        raise ValueError("Invalid scPCA input dimensions or condition labels")
    n, p = x_np.shape
    k = min(args.rank, n - 1, p - 1)
    q = len(conditions) - 1
    lib_np = np.maximum(x_np.sum(axis=1), 1.0)
    loglib_np = np.log(lib_np)

    x = torch.from_numpy(x_np)
    phi = torch.from_numpy(np.clip(phi_np, 1e-6, 100.0))
    loglib = torch.from_numpy(loglib_np)
    cond_t = torch.from_numpy(cond.astype(np.int64))


    design = np.zeros((len(conditions), q + 1), dtype=np.float32)
    design[:, 0] = 1.0
    for j in range(q):
        design[j + 1, j + 1] = 1.0
    design /= np.sqrt((design**2).sum(axis=1, keepdims=True))
    design_t = torch.from_numpy(design)

    gene_prop = (x_np.sum(axis=0) + 0.1) / (lib_np.sum() + 0.1 * p)
    base_intercept = np.log(np.maximum(gene_prop, 1e-10)).astype(np.float32)
    z = torch.nn.Parameter(0.05 * torch.randn(n, k))
    w = torch.nn.Parameter(0.05 * torch.randn(q + 1, k, p))
    v_init = np.tile(base_intercept, (len(conditions), 1))
    v = torch.nn.Parameter(torch.from_numpy(v_init))
    optimizer = torch.optim.Adam([z, w, v], lr=0.01)
    batch_size = min(args.batch_size, n)
    losses = []
    stable = 0
    for step in range(args.steps):
        idx = torch.randperm(n)[:batch_size]
        state_w = torch.einsum("bd,dkp->bkp", design_t[cond_t[idx]], w)
        log_mu = loglib[idx, None] + v[cond_t[idx]] + torch.einsum(
            "bk,bkp->bp", z[idx], state_w
        )
        mu = torch.exp(torch.clamp(log_mu, max=15.0))
        likelihood = nb_nll(x[idx], mu, phi[None, :]).mean() * n
        prior = 0.5 * (w.square().sum() + v.square().sum()) / p
        prior = prior + 0.5 * z.square().sum() / (0.1**2 * n)
        loss = likelihood + prior
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_([z, w, v], 100.0)
        optimizer.step()
        value = float(loss.detach())
        losses.append(value)
        if step > 50:
            old = np.mean(losses[-41:-21])
            new = np.mean(losses[-20:])
            rel = abs(old - new) / max(1.0, abs(old))
            stable = stable + 1 if rel < 1e-5 else 0
            if stable >= 20:
                break

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mu0 = torch.from_numpy(mu0_np)
    denom = torch.sqrt(mu0 + mu0.square() * phi[None, :]).clamp_min(1e-8)
    with torch.no_grad():
        for j in range(q):
            wr = torch.einsum("d,dkp->kp", design_t[0], w)
            wt = torch.einsum("d,dkp->kp", design_t[j + 1], w)
            log_ref = loglib[:, None] + v[0] + z @ wr
            log_tgt = loglib[:, None] + v[j + 1] + z @ wt
            mu_ref = torch.exp(torch.clamp(log_ref, max=15.0))
            mu_tgt = torch.exp(torch.clamp(log_tgt, max=15.0))
            r_ref = torch.clamp((mu_ref - mu0) / denom, -args.clip, args.clip)
            r_tgt = torch.clamp((mu_tgt - mu0) / denom, -args.clip, args.clip)
            np.save(
                outdir / f"response_{j + 1}.npy",
                (r_tgt - r_ref).cpu().numpy().astype(np.float64),
            )
    with open(outdir / "diagnostics.json", "w", encoding="utf-8") as handle:
        json.dump(
            {
                "steps": len(losses),
                "converged": stable >= 20,
                "loss_initial": losses[0],
                "loss_final": losses[-1],
                "rank": k,
            },
            handle,
            indent=2,
        )


if __name__ == "__main__":
    main()
