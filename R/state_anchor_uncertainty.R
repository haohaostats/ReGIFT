# Frozen r4 conditional state-deviation shrinkage, integrated in version 0.1.0.
# Conditional on the fixed working response, donor jackknife estimates the
# variance of state deviations from the common response anchor.
# Working model: d_hat | d ~ N(d,v), d ~ N(0,tau2), with shared tau2 estimated
# by positive-part method of moments over supported state/gene deviations.
# This is a point-estimation rule, not calibrated inferential intervals.
.state_deviation_posterior <- function(delta, variance) {
  valid <- is.finite(delta) & is.finite(variance) & variance >= 0
  if (!any(valid)) return(list(mean=delta,weight=matrix(NA_real_,nrow(delta),ncol(delta)),
    tau2=NA_real_,status='no finite variance estimates'))
  tau2 <- max(0,mean(delta[valid]^2-variance[valid]))
  weight <- matrix(0,nrow(delta),ncol(delta),dimnames=dimnames(delta))
  if (tau2 > 0) weight[valid] <- tau2/(tau2+variance[valid])
  # Infinite variance is an explicit global-response borrowing rule.
  posterior <- delta*weight
  posterior[!is.finite(posterior)] <- 0
  list(mean=posterior,weight=weight,tau2=tau2,status='conditional normal-means moment estimate')
}
.regift_state_uncertainty <- function(Y,meta,contrasts,A0,A_state,ridge,
                                      condition_main_strata=NULL,condition_main=TRUE) {
  base<-asNamespace('ReGIFT')
  unchanged<-function(status)list(A_state=A_state,status=status,applied=FALSE)
  if(is.null(A_state)||!isTRUE(condition_main))return(unchanged('no active condition/state anchors'))
  if(nrow(contrasts)!=1L)return(unchanged('state shrinkage supports one contrast only'))
  donors<-unique(as.character(meta$donor));states<-dimnames(A_state)[[3]]
  if(is.null(states))states<-unique(as.character(meta$state))
  if(length(donors)<3L)return(unchanged('insufficient donors for conditional jackknife'))
  support<-attr(A_state,'state_anchor_support')
  if(is.null(support))stop('State support metadata is required by the state shrinkage.')
  p<-ncol(Y);D<-length(donors);H<-length(states)
  deleted<-array(NA_real_,c(D,p,H));rank<-matrix(NA_integer_,D,H)
  design_by_donor<-vapply(donors,function(s)length(unique(meta$condition[meta$donor==s])),integer(1))
  paired<-all(design_by_donor>1L)
  if(!paired&&any(design_by_donor>1L))return(unchanged('state shrinkage requires a purely paired or unpaired design'))
  jack_groups<-if(paired)list(paired=seq_along(donors))else
    split(seq_along(donors),vapply(donors,function(s)as.character(meta$condition[which(meta$donor==s)[1L]]),character(1)))
  failures<-character()
  for(s in seq_along(donors)) {
    ii<-which(meta$donor!=donors[s]);m<-meta[ii,,drop=FALSE]
    available<-states[states%in%m$state]
    ans<-tryCatch({
      strata<-condition_main_strata
      if(!is.null(strata)&&length(strata)&&!is.character(strata))strata<-strata[ii]
      w<-base$.regift_condition_main_weights(m,strata)
      d<-base$.regift_cell_design(m,contrasts)
      a<-t(base$.regift_condition_slopes(d,Y[ii,,drop=FALSE],w,ridge,
        strata=base$.regift_anchor_baseline(m,TRUE)))
      ast<-base$.regift_robust_state_effects(Y[ii,,drop=FALSE],m,d,a,available,ridge)
      list(ast=ast,support=attr(ast,'state_anchor_support'))
    },error=function(e)e)
    if(inherits(ans,'error')) {
      failures<-c(failures,paste(donors[s],conditionMessage(ans),sep=': '));next
    }
    for(h in seq_along(available)) {
      hh<-match(available[h],states)
      deleted[s,,hh]<-ans$ast[,1,h]
      rank[s,hh]<-ans$support$estimable_rank[h]
    }
  }
  if(length(failures))return(c(unchanged('global jackknife fit unavailable'),list(failures=failures)))
  delta<-matrix(A_state[,1,],p,H,dimnames=list(colnames(Y),states))
  variance<-matrix(Inf,p,H,dimnames=dimnames(delta))
  stable<-logical(H)
  for(h in seq_len(H)) {
    stable[h]<-support$estimable_rank[h]==1L && all(rank[,h]==1L,na.rm=FALSE) && all(is.finite(deleted[,,h]))
    if(isTRUE(stable[h])) {
      variance[,h]<-0
      for(jj in jack_groups) {
        z<-matrix(deleted[jj,,h,drop=FALSE],length(jj),p)
        variance[,h]<-variance[,h]+(length(jj)-1)/length(jj)*
          colSums(sweep(z,2,colMeans(z),'-')^2)
      }
    }
  }
  posterior<-.state_deviation_posterior(delta,variance)
  if(!is.finite(posterior$tau2))return(unchanged('no stable state support for prior variance'))
  out<-A_state;out[,1,]<-posterior$mean
  support$uncertainty_rule<-ifelse(stable,'conditional donor-jackknife shrinkage',
    'support changes on donor deletion; deviation pooled to common anchor')
  attr(out,'state_anchor_support')<-support
  list(A_state=out,status=posterior$status,applied=TRUE,tau2=posterior$tau2,
    raw_delta=delta,conditional_variance=variance,weight=posterior$weight,
    stable_support=setNames(stable,states),deleted_rank=rank,donors=donors,
    jackknife_design=if(paired) "paired donor deletion" else "condition-stratified donor deletion",
    estimand_note='conditional-on-working-response variance; robust locations and support changes limit jackknife regularity')
}
