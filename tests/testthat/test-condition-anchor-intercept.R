test_that('condition anchors separate unequal-arm state baselines from effects', {
  meta<-expand.grid(donor=paste0('D',1:8),state=c('S1','S2'),cell=1:6,stringsAsFactors=FALSE)
  meta$condition<-ifelse(meta$donor%in%paste0('D',1:5),'ctrl','stim')
  meta$sample<-paste(meta$donor,meta$condition,sep='_')
  contrasts<-regift_contrasts(c('ctrl','stim'))
  d<-ReGIFT:::.regift_cell_design(meta,contrasts)
  genes<-seq(.5,1.5,length.out=8)
  for(effect in c(0,.8)) {
    Y<-(ifelse(meta$state=='S1',2,-2)+d[,1]*effect)%o%genes
    w<-ReGIFT:::.regift_condition_main_weights(meta)
    A<-t(ReGIFT:::.regift_condition_slopes(d,Y,w,ridge=1e-10))
    Ast<-ReGIFT:::.regift_robust_state_effects(Y,meta,d,A,c('S1','S2'),ridge=1e-10)
    expect_equal(A[,1]+Ast[,1,1],effect*genes,tolerance=1e-8)
    expect_equal(A[,1]+Ast[,1,2],effect*genes,tolerance=1e-8)
  }
})

test_that('condition slopes agree with weighted regression including an intercept', {
  x<-cbind(seq(-1,1,length.out=30),rep(c(-1,1),15))
  y<-cbind(sin(1:30),cos(1:30));w<-seq(.1,2,length.out=30)
  reference<-qr.solve(cbind(1,x)*sqrt(w),y*sqrt(w))[-1,,drop=FALSE]
  observed<-ReGIFT:::.regift_condition_slopes(x,y,w,ridge=0)
  expect_equal(unname(observed),unname(reference),tolerance=1e-10)
  expect_equal(ReGIFT:::.regift_condition_slopes(x,sweep(y,2,c(12,-4),'+'),w,0),observed,tolerance=1e-10)
  expect_error(ReGIFT:::.regift_condition_slopes(matrix(1,30,1),y,w),'not estimable')
})
test_that('unsupported state directions borrow global response with recorded support', {
  x<-cbind(rep(c(-1,1),10),0)
  y<-cbind(3+2*x[,1],-4-x[,1]);w<-rep(1/20,20)
  info<-ReGIFT:::.regift_condition_slopes(x,y,w,0,TRUE,TRUE)
  expect_equal(info$rank,1L)
  expect_true(info$partial_pooling)
  expect_equal(unname(info$coefficients),rbind(c(2,-1),c(0,0)),tolerance=1e-10)
  expect_equal(unname(info$identified_projection),diag(c(1,0)),tolerance=1e-10)
  single<-ReGIFT:::.regift_condition_slopes(matrix(.5,20,1),y,w,0,TRUE,TRUE)
  expect_equal(single$rank,0L)
  expect_equal(single$coefficients,matrix(0,1,2))
})

test_that('pooled condition anchors are invariant to unsupported state baseline shifts', {
  meta<-expand.grid(donor=paste0('D',1:6),cell=1:6,stringsAsFactors=FALSE)
  meta$condition<-ifelse(meta$donor%in%paste0('D',1:3),'ctrl','stim')
  meta$sample<-meta$donor;meta$state<-'shared'
  extra<-meta[meta$condition=='ctrl',];extra$state<-'rare'
  meta<-rbind(meta,extra)
  cc<-regift_contrasts(c('ctrl','stim'));d<-ReGIFT:::.regift_cell_design(meta,cc)
  g<-seq(.5,1.5,length.out=8)
  for(shift in c(0,20)) {
    y<-(1+d[,1]*.8+ifelse(meta$state=='rare',shift,0))%o%g
    fit<-regift_fit(y,meta,cc,K=2,H=0,lambda_B=.1,max_iter=1,
      z_backend='blocked',svd_backend='dense')
    expect_equal(as.numeric(fit$A),.8*g,tolerance=1e-4)
    support<-fit$state_anchor_support
    rare<-which(support$state=='rare')
    expect_true(support$partial_pooling[rare])
    expect_equal(support$estimable_rank[rare],0L)
    h<-match('rare',fit$state_levels)
    expect_equal(as.numeric(fit$A_state[,,h]),rep(0,8),tolerance=1e-10)
    expect_equal(support$identified_projection[[rare]],matrix(0,1,1))
  }
  confounded<-ifelse(meta$condition=='ctrl','only_ctrl','only_stim')
  expect_error(ReGIFT:::.regift_condition_slopes(d,y,rep(1/nrow(meta),nrow(meta)),
    strata=confounded),'not estimable')
})