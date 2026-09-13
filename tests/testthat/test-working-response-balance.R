test_that("working-response exposure reference is donor balanced", {
  counts <- matrix(c(10, 20, 30, 40,
                     12, 22, 32, 42,
                     50, 60, 70, 80,
                     55, 65, 75, 85), nrow = 4, byrow = TRUE)
  meta <- data.frame(donor = c("D1", "D1", "D2", "D2"),
                     sample = paste0("S", 1:4))
  base <- ReGIFT::regift_working_response(counts, meta)
  duplicated <- ReGIFT::regift_working_response(rbind(counts, counts[1:2, ]),
    rbind(meta, meta[1:2, ]))
  expect_equal(base$library_geomean, duplicated$library_geomean,
               tolerance = 1e-12)
  expect_equal(base$rate, duplicated$rate, tolerance = 1e-12)
})

test_that("working representation is invariant to unequal-donor row order", {
  set.seed(9274001)
  counts <- matrix(rpois(120,10),6,20)
  meta <- data.frame(donor=c("B","B","A","A","A","A"),
    sample=c("B0","B1","A0","A0","A1","A1"),
    condition=c("C0","C1","C0","C0","C1","C1"))
  a <- regift_working_response(counts,meta)
  order <- c(3:6,1:2)
  b <- regift_working_response(counts[order,],meta[order,])
  expect_equal(a$library_geomean,b$library_geomean,tolerance=1e-12)
  expect_equal(a$rate,b$rate,tolerance=1e-12)
  expect_equal(a$Y[order,],b$Y,tolerance=1e-12)
  expect_equal(a$phi,b$phi,tolerance=1e-12)
  expected_geomean <- exp(mean(c(mean(log(rowSums(counts)[1:2])),
    mean(log(rowSums(counts)[3:6])))))
  expect_equal(a$library_geomean,expected_geomean,tolerance=1e-12)
})