test_that("constrained rescue preserves multi-condition geometry and penalties", {
  for (h in c(0L,1L)) {
    sim <- regift_simulate(A=3,S=4,p=24,cells=10,K0=2,H=h,seed=772+h)
    wr <- regift_working_response(sim$counts,sim$meta)
    f <- regift_fit(wr$Y,sim$meta,sim$contrasts,K=2,H=h,max_iter=3,
      lambda_B=.001,lambda_Delta=.1)
    d <- ReGIFT:::.regift_cell_design(f$meta,f$contrasts)
    obj <- function(x)ReGIFT:::.regift_objective(x,wr$Y,d,
      f$lambda_B,f$lambda_Delta,f$sigma,1.345)
    ans <- ReGIFT:::.regift_constrained_rescue(f,wr$Y,d,f$sigma,1.345,obj,
      f$lambda_B,f$lambda_Delta,256L)
    z <- ans$fit
    expect_lte(ans$objective,obj(f)+1e-12)
    nw <- z$normalization_weights; nw <- nw/sum(nw)
    expect_equal(crossprod(z$Z,z$Z*nw),diag(2),tolerance=1e-7)
    for(q in seq_len(ncol(d))) {
      inf <- vapply(z$donors,function(s)any(abs(d[z$meta$donor==s,q])>0),logical(1))
      expect_lt(max(abs(apply(z$Delta[inf,,,q,drop=FALSE],c(2,3),mean))),1e-7)
    }
    for(s in seq_along(z$donors)) if(h) {
      ii <- which(z$meta$donor==z$donors[s]); zz <- z$Z[ii,,drop=FALSE]
      bio <- cbind(d[ii,,drop=FALSE],zz,
        do.call(cbind,lapply(seq_len(ncol(d)),function(q)zz*d[ii,q])))
      expect_lt(max(abs(crossprod(bio,z$U[[s]]))),1e-7)
    }
    if(ans$block=='Z_tangent') {
      expect_identical(z$B,f$B)
      expect_identical(z$Delta,f$Delta)
      expect_identical(z$L0,f$L0)
    }
  }
})
