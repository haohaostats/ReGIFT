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
