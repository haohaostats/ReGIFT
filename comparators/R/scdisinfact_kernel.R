scdisinfact_kernel_fit <- function(counts, meta, contrasts, working,
                                    shared_rank=10L, condition_rank=5L, hidden=128L,
                                    steps=1000L, batch_size=1024L, seed=20260829L,
                                    python=Sys.getenv("REGIFT_COMPARATOR_PYTHON", "python"),
                                    script=.comparator_script("scdisinfact_kernel.py")) {
  stopifnot(requireNamespace("RcppCNPy",quietly=TRUE),file.exists(script))
  work<-tempfile("scdisinfact_",tmpdir=tempdir()); dir.create(work,recursive=TRUE)
  on.exit(unlink(work,recursive=TRUE,force=TRUE),add=TRUE)
  f<-file.path(work,c("counts.npy","mu0.npy","phi.npy","meta.csv"))
  RcppCNPy::npySave(f[1],as.matrix(counts)); RcppCNPy::npySave(f[2],as.matrix(working$mu0)); RcppCNPy::npySave(f[3],as.numeric(working$phi)); write.csv(meta,f[4],row.names=FALSE)
  outdir<-file.path(work,"out")
  args<-c(script,"--counts",f[1],"--meta",f[4],"--mu0",f[2],"--phi",f[3],"--clip",working$clip,"--conditions",paste(colnames(contrasts),collapse=","),"--outdir",outdir,"--shared-rank",shared_rank,"--condition-rank",condition_rank,"--hidden",hidden,"--steps",steps,"--batch-size",batch_size,"--seed",seed)
  status<-system2(python,args = shQuote(as.character(args))); if(!identical(status,0L)) stop("scDisInFact kernel failed with status ",status)
  response<-lapply(seq_len(nrow(contrasts)),function(q) RcppCNPy::npyLoad(file.path(outdir,paste0("response_",q,".npy"))))
  names(response)<-rownames(contrasts)
  list(population_response=response,shared_latent=RcppCNPy::npyLoad(file.path(outdir,"shared_latent.npy")),diagnostics=jsonlite::fromJSON(file.path(outdir,"diagnostics.json")),extraction="in-house scDisInFact disentangled NB-VAE core with condition-prior counterfactual decoding")
}
