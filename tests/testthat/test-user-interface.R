test_that("regift high-level interface fits counts and returns summaries", {
  sim <- regift_simulate(A=2,S=4,p=30,cells=10,K0=2,H=0,seed=1601)
  sim$meta$condition <- ifelse(sim$meta$condition=="C0","control","treated")
  ans <- regift(sim$counts,sim$meta,state="state",reference="control",
                K=2,H=0,max_iter=6,tol=1,threads=1)
  expect_s3_class(ans,"regift_analysis")
  expect_s3_class(ans$fit,"regift_fit")
  expect_named(ans$population_response,"treated-control")
  expect_equal(dim(ans$population_response[[1]]),dim(sim$counts))
  tab <- regift_response_table(ans)
  expect_true(all(c("contrast","state","gene","mean_response",
                    "rms_response","direction") %in% names(tab)))
  expect_equal(nrow(summary(ans)),1L)
})

test_that("regift validates sample mappings and reference condition", {
  sim <- regift_simulate(A=2,S=3,p=20,cells=10,K0=2,H=0,seed=1602)
  bad <- sim$meta
  bad$donor[bad$sample==bad$sample[1]][1] <- "wrong"
  expect_error(regift(sim$counts,bad,K=2,H=0,max_iter=2),
               "Each sample must map")
  expect_error(regift(sim$counts,sim$meta,reference="missing",K=2,H=0,max_iter=2),
               "reference")
})

test_that("bundled example has the documented user contract", {
  data("regift_example",package="ReGIFT")
  expect_s4_class(regift_example$counts,"dgCMatrix")
  expect_equal(nrow(regift_example$counts),nrow(regift_example$meta))
  expect_true(all(c("cell","donor","sample","condition","state") %in%
                    names(regift_example$meta)))
  expect_setequal(unique(regift_example$meta$condition),c("control","T1D"))
  expect_equal(length(unique(regift_example$meta$donor)),10L)
  expect_equal(ncol(regift_example$counts),240L)
  expect_equal(regift_example$provenance$accession,"GSE148073")
})
