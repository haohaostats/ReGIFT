test_that('rare paired-state anchors block donor baselines', {
  meta<-expand.grid(donor=paste0('D',1:4),condition=c('ctrl','stim'),cell=1:4,stringsAsFactors=FALSE)
  meta$sample<-paste(meta$donor,meta$condition,sep='_');meta$state<-'common'
  rare<-meta[!(meta$donor=='D3'&meta$condition=='stim') & !(meta$donor=='D4'&meta$condition=='ctrl'),]
  rare$state<-'rare';meta<-rbind(meta,rare)
  cc<-regift_contrasts(c('ctrl','stim'));d<-ReGIFT:::.regift_cell_design(meta,cc)
  beta<-seq(.4,1.2,length.out=8)
  for(baseline in list(c(D1=2,D2=-3,D3=10,D4=-10),c(D1=13,D2=-10,D3=110,D4=-70))) {
    y<-(5+baseline[meta$donor])%o%rep(1,8)+d%*%t(beta)
    fit<-regift_fit(y,meta,cc,K=2,H=0,lambda_B=.1,max_iter=1,svd_backend='dense')
    for(h in seq_along(fit$state_levels))
      expect_equal(as.numeric(fit$A+fit$A_state[,,h]),beta,tolerance=1e-4)
    expect_identical(fit$anchor_baseline_model,'within-donor')
    expect_true(all(fit$state_anchor_support$baseline_model=='within-donor'))
  }
})

test_that('connected partial pairs identify multiple contrasts within donors', {
  base<-data.frame(donor=rep(paste0('D',1:3),each=2),
    condition=c('C0','C1','C1','C2','C0','C2'),stringsAsFactors=FALSE)
  meta<-base[rep(seq_len(nrow(base)),each=8),]
  meta$sample<-paste(meta$donor,meta$condition,sep='_')
  meta$state<-rep(rep(c('s1','s2'),each=4),nrow(base))
  cc<-regift_contrasts(c('C0','C1','C2'));d<-ReGIFT:::.regift_cell_design(meta,cc)
  beta<-rbind(seq(.1,.8,length.out=8),seq(-.9,-.2,length.out=8))
  baseline<-c(D1=5,D2=-10,D3=20)[meta$donor]+ifelse(meta$state=='s1',2,-2)
  y<-baseline%o%rep(1,8)+d%*%beta
  fit<-regift_fit(y,meta,cc,K=2,H=0,lambda_B=.1,max_iter=1,svd_backend='dense')
  for(h in seq_along(fit$state_levels))
    expect_equal(unname(fit$A+fit$A_state[,,h]),t(beta),tolerance=1e-4)
  expect_true(all(fit$state_anchor_support$paired_donors==0L))
  expect_true(all(fit$state_anchor_support$estimable_rank==2L))
})