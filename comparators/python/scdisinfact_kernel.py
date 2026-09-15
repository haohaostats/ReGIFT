import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn

def nb_nll(x, mu, phi):
    r = 1.0 / phi
    return -(torch.lgamma(x+r)-torch.lgamma(r)-torch.lgamma(x+1)+
             r*(torch.log(r)-torch.log(r+mu))+
             x*(torch.log(mu.clamp_min(1e-8))-torch.log(r+mu)))

class Model(nn.Module):
    def __init__(self,p,nb,nc,ks,ku,h):
        super().__init__()
        self.enc=nn.Sequential(nn.Linear(p,h),nn.Softplus(),nn.Linear(h,h),nn.Softplus())
        self.sm=nn.Linear(h,ks); self.sv=nn.Linear(h,ks)
        self.um=nn.Linear(h,ku); self.uv=nn.Linear(h,ku)
        self.prior_mean=nn.Parameter(torch.randn(nc,ku)*.05)
        self.prior_logvar=nn.Parameter(torch.zeros(nc,ku))
        self.classifier=nn.Linear(ku,nc)
        self.batch=nn.Embedding(nb,ku)
        self.dec=nn.Sequential(nn.Linear(ks+2*ku,h),nn.Softplus(),nn.Linear(h,p))
    def encode(self,x):
        h=self.enc(x)
        return self.sm(h),self.sv(h).clamp(-8,8),self.um(h),self.uv(h).clamp(-8,8)
    def decode(self,zs,zu,b,lib):
        logits=self.dec(torch.cat([zs,zu,self.batch(b)],1))
        return lib[:,None]*torch.softmax(logits,1)

def main():
    a=argparse.ArgumentParser()
    for k in ["counts","meta","mu0","phi","conditions","outdir"]: a.add_argument("--"+k,required=True)
    a.add_argument("--clip",type=float,required=True); a.add_argument("--shared-rank",type=int,default=10)
    a.add_argument("--condition-rank",type=int,default=5); a.add_argument("--hidden",type=int,default=128)
    a.add_argument("--steps",type=int,default=1000); a.add_argument("--batch-size",type=int,default=1024)
    a.add_argument("--seed",type=int,default=20260829); args=a.parse_args()
    torch.manual_seed(args.seed); np.random.seed(args.seed)
    xn=np.load(args.counts).astype("float32"); mu0n=np.load(args.mu0).astype("float32")
    phin=np.load(args.phi).astype("float32").reshape(-1); meta=pd.read_csv(args.meta,dtype=str)
    conditions=args.conditions.split(","); cc=pd.Categorical(meta.condition,categories=conditions).codes
    bc,bnames=pd.factorize(meta["sample"],sort=True); n,p=xn.shape
    libn=np.maximum(xn.sum(1),1).astype("float32"); en=np.log1p(1e4*xn/libn[:,None]).astype("float32")
    x=torch.from_numpy(xn); enc=torch.from_numpy(en); lib=torch.from_numpy(libn)
    cond=torch.from_numpy(cc.astype("int64")); batch=torch.from_numpy(bc.astype("int64"))
    phi=torch.from_numpy(np.clip(phin,1e-6,100)); model=Model(p,len(bnames),len(conditions),args.shared_rank,args.condition_rank,args.hidden)
    opt=torch.optim.Adam(model.parameters(),lr=1e-3); losses=[]; stable=0; bs=min(args.batch_size,n)
    for step in range(args.steps):
        ii=torch.randperm(n)[:bs]; sm,sv,um,uv=model.encode(enc[ii])
        zs=sm+torch.randn_like(sm)*torch.exp(.5*sv); zu=um+torch.randn_like(um)*torch.exp(.5*uv)
        mu=model.decode(zs,zu,batch[ii],lib[ii]); recon=nb_nll(x[ii],mu,phi[None]).sum(1).mean()
        kls=-.5*(1+sv-sm.square()-sv.exp()).sum(1).mean()
        pm=model.prior_mean[cond[ii]]; pv=model.prior_logvar[cond[ii]].clamp(-8,8)
        klu=.5*(pv-uv+(uv.exp()+(um-pm).square())/pv.exp()-1).sum(1).mean()
        cls=nn.functional.cross_entropy(model.classifier(um),cond[ii])
        loss=recon+kls+klu+cls
        opt.zero_grad(); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),100); opt.step()
        losses.append(float(loss.detach()))
        if step>100:
            old,new=np.mean(losses[-61:-31]),np.mean(losses[-30:]); stable=stable+1 if abs(old-new)/max(1,abs(old))<1e-4 else 0
            if stable>=20: break
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    with torch.no_grad():
        zs,_,_,_=model.encode(enc); mu0=torch.from_numpy(mu0n); den=torch.sqrt(mu0+mu0.square()*phi[None]).clamp_min(1e-8)
        decoded=[]
        for c in range(len(conditions)):
            acc=torch.zeros((n,p)); zu=model.prior_mean[c][None].expand(n,-1)
            for b in range(len(bnames)):
                bb=torch.full((n,),b,dtype=torch.int64); acc+=model.decode(zs,zu,bb,lib)
            decoded.append(acc/len(bnames))
        for q in range(1,len(conditions)):
            rr=torch.clamp((decoded[q]-mu0)/den,-args.clip,args.clip)-torch.clamp((decoded[0]-mu0)/den,-args.clip,args.clip)
            np.save(out/f"response_{q}.npy",rr.numpy().astype("float64"))
        np.save(out/"shared_latent.npy",zs.numpy().astype("float64"))
    with open(out/"diagnostics.json","w") as f: json.dump({"steps":len(losses),"converged":stable>=20,"loss_initial":losses[0],"loss_final":losses[-1]},f)
if __name__=="__main__": main()
