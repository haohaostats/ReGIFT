import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn


def nb_nll(x, mu, phi):
    size = 1.0 / phi
    return -(
        torch.lgamma(x + size) - torch.lgamma(size) - torch.lgamma(x + 1.0)
        + size * (torch.log(size) - torch.log(size + mu))
        + x * (torch.log(mu.clamp_min(1e-8)) - torch.log(size + mu))
    )


class ScVICore(nn.Module):
    def __init__(self, genes, batches, latent, hidden):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(genes, hidden), nn.Softplus(),
                                     nn.Linear(hidden, hidden), nn.Softplus())
        self.z_mean = nn.Linear(hidden, latent)
        self.z_logvar = nn.Linear(hidden, latent)
        self.batch = nn.Embedding(batches, latent)
        self.decoder = nn.Sequential(nn.Linear(2 * latent, hidden), nn.Softplus(),
                                     nn.Linear(hidden, genes))

    def encode(self, x):
        h = self.encoder(x)
        return self.z_mean(h), self.z_logvar(h).clamp(-8.0, 8.0)

    def decode(self, z, batch, library):
        logits = self.decoder(torch.cat([z, self.batch(batch)], dim=1))
        return library[:, None] * torch.softmax(logits, dim=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--mu0", required=True)
    ap.add_argument("--phi", required=True)
    ap.add_argument("--clip", required=True, type=float)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rank", type=int, default=10)
    ap.add_argument("--hidden", type=int, default=128)
    ap.add_argument("--steps", type=int, default=1000)
    ap.add_argument("--batch-size", type=int, default=1024)
    ap.add_argument("--seed", type=int, default=20260829)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    x_np = np.load(args.counts).astype(np.float32)
    mu0_np = np.load(args.mu0).astype(np.float32)
    phi_np = np.load(args.phi).astype(np.float32).reshape(-1)
    meta = pd.read_csv(args.meta, dtype=str)
    batch_codes, batch_names = pd.factorize(meta["sample"], sort=True)
    n, p = x_np.shape
    if mu0_np.shape != x_np.shape or len(meta) != n:
        raise ValueError("Invalid scVI input dimensions")
    lib_np = np.maximum(x_np.sum(1), 1.0).astype(np.float32)
    enc_np = np.log1p(1e4 * x_np / lib_np[:, None]).astype(np.float32)
    x = torch.from_numpy(x_np)
    enc = torch.from_numpy(enc_np)
    lib = torch.from_numpy(lib_np)
    batch = torch.from_numpy(batch_codes.astype(np.int64))
    phi = torch.from_numpy(np.clip(phi_np, 1e-6, 100.0))
    model = ScVICore(p, len(batch_names), min(args.rank, p - 1), args.hidden)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    losses, stable = [], 0
    bs = min(args.batch_size, n)
    for step in range(args.steps):
        idx = torch.randperm(n)[:bs]
        zm, zv = model.encode(enc[idx])
        z = zm + torch.randn_like(zm) * torch.exp(0.5 * zv)
        mu = model.decode(z, batch[idx], lib[idx])
        recon = nb_nll(x[idx], mu, phi[None, :]).sum(1).mean()
        kl = -0.5 * (1 + zv - zm.square() - zv.exp()).sum(1).mean()
        loss = recon + kl
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 100.0)
        opt.step()
        losses.append(float(loss.detach()))
        if step > 100:
            old, new = np.mean(losses[-61:-31]), np.mean(losses[-30:])
            stable = stable + 1 if abs(old - new) / max(1.0, abs(old)) < 1e-4 else 0
            if stable >= 20:
                break

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    with torch.no_grad():
        z, _ = model.encode(enc)


        rho = torch.zeros((n, p))
        unit_library = torch.ones(n)
        for b in range(len(batch_names)):
            bb = torch.full((n,), b, dtype=torch.int64)
            rho += model.decode(z, bb, unit_library)
        rho /= len(batch_names)
        adjusted_mu = lib[:, None] * rho
        mu0 = torch.from_numpy(mu0_np)
        denom = torch.sqrt(mu0 + mu0.square() * phi[None, :]).clamp_min(1e-8)
        adjusted = torch.clamp((adjusted_mu - mu0) / denom,
                               -args.clip, args.clip)
        np.save(outdir / "adjusted.npy", adjusted.numpy().astype(np.float64))
        np.save(outdir / "latent.npy", z.numpy().astype(np.float64))
    with open(outdir / "diagnostics.json", "w", encoding="utf-8") as f:
        json.dump({"steps": len(losses), "converged": stable >= 20,
                   "loss_initial": losses[0], "loss_final": losses[-1],
                   "rank": model.z_mean.out_features,
                   "batches": len(batch_names)}, f, indent=2)


if __name__ == "__main__":
    main()
