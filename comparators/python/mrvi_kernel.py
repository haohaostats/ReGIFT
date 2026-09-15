import argparse,json
from pathlib import Path
import numpy as np,pandas as pd,torch
from torch import nn
def nb(x,mu,phi):
 r=1/phi; return -(torch.lgamma(x+r)-torch.lgamma(r)-torch.lgamma(x+1)+r*(torch.log(r)-torch.log(r+mu))+x*(torch.log(mu.clamp_min(1e-8))-torch.log(r+mu)))
class M(nn.Module):
 def __init__(self,p,ns,k,h):
  super().__init__(); self.e=nn.Sequential(nn.Linear(p,h),nn.Softplus(),nn.Linear(h,h),nn.Softplus()); self.um=nn.Linear(h,k); self.uv=nn.Linear(h,k); self.se=nn.Embedding(ns,k); self.shift=nn.Sequential(nn.Linear(2*k,h),nn.ReLU(),nn.Linear(h,k)); self.d=nn.Sequential(nn.Linear(k,h),nn.Softplus(),nn.Linear(h,p))
 def enc(self,x): h=self.e(x); return self.um(h),self.uv(h).clamp(-8,8)
 def z(self,u,s): return u+self.shift(torch.cat([u,self.se(s)],1))
 def dec(self,z,l): return l[:,None]*torch.softmax(self.d(z),1)
def main():
 a=argparse.ArgumentParser()
 for k in ["counts","meta","mu0","phi","conditions","outdir"]: a.add_argument("--"+k,required=True)
 a.add_argument("--clip",type=float,required=True); a.add_argument("--rank",type=int,default=10); a.add_argument("--hidden",type=int,default=128); a.add_argument("--steps",type=int,default=1000); a.add_argument("--batch-size",type=int,default=1024); a.add_argument("--seed",type=int,default=20260829); x=a.parse_args()
 torch.manual_seed(x.seed); np.random.seed(x.seed); Xn=np.load(x.counts).astype("float32"); M0n=np.load(x.mu0).astype("float32"); phin=np.load(x.phi).astype("float32").reshape(-1); md=pd.read_csv(x.meta,dtype=str); conds=x.conditions.split(","); sc,snames=pd.factorize(md["sample"],sort=True); sample_cond=md.drop_duplicates("sample").set_index("sample").loc[snames,"condition"].values; n,p=Xn.shape; ln=np.maximum(Xn.sum(1),1).astype("float32"); en=np.log1p(1e4*Xn/ln[:,None]).astype("float32")
 X=torch.from_numpy(Xn); E=torch.from_numpy(en); L=torch.from_numpy(ln); S=torch.from_numpy(sc.astype("int64")); phi=torch.from_numpy(np.clip(phin,1e-6,100)); model=M(p,len(snames),min(x.rank,p-1),x.hidden); opt=torch.optim.Adam(model.parameters(),lr=1e-3); losses=[]; stable=0; bs=min(x.batch_size,n)
 for step in range(x.steps):
  ii=torch.randperm(n)[:bs]; um,uv=model.enc(E[ii]); u=um+torch.randn_like(um)*torch.exp(.5*uv); z=model.z(u,S[ii]); mu=model.dec(z,L[ii]); recon=nb(X[ii],mu,phi[None]).sum(1).mean(); kl=-.5*(1+uv-um.square()-uv.exp()).sum(1).mean(); proximity=.5*(z-u).square().sum(1).mean(); loss=recon+kl+proximity; opt.zero_grad(); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),100); opt.step(); losses.append(float(loss.detach()))
  if step>100:
   old,new=np.mean(losses[-61:-31]),np.mean(losses[-30:]); stable=stable+1 if abs(old-new)/max(1,abs(old))<1e-4 else 0
   if stable>=20: break
 out=Path(x.outdir); out.mkdir(parents=True,exist_ok=True)
 with torch.no_grad():
  u,_=model.enc(E); mu0=torch.from_numpy(M0n); den=torch.sqrt(mu0+mu0.square()*phi[None]).clamp_min(1e-8); decoded=[]
  for c in conds:
   ids=np.where(sample_cond==c)[0]; acc=torch.zeros((n,p))
   for s in ids: acc+=model.dec(model.z(u,torch.full((n,),int(s),dtype=torch.int64)),L)
   decoded.append(acc/max(len(ids),1))
  for q in range(1,len(conds)):
   rr=torch.clamp((decoded[q]-mu0)/den,-x.clip,x.clip)-torch.clamp((decoded[0]-mu0)/den,-x.clip,x.clip); np.save(out/f"response_{q}.npy",rr.numpy().astype("float64"))
  np.save(out/"u.npy",u.numpy().astype("float64"))
 with open(out/"diagnostics.json","w") as f: json.dump({"steps":len(losses),"converged":stable>=20,"loss_initial":losses[0],"loss_final":losses[-1]},f)
if __name__=="__main__": main()
