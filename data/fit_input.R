# Working representation and fitting settings for the published input matrices.
regift_input_working <- function(x, dataset) {
  X <- as.matrix(x$counts)
  colnames(X) <- as.character(x$genes)
  if (dataset != 'parse_tnf') return(ReGIFT::regift_working_response(X, x$meta)$Y)
  m <- x$meta
  stopifnot(length(m$library_size) == nrow(X), all(m$library_size > 0))
  groups <- split(seq_len(nrow(X)), factor(m$donor))
  reference <- exp(mean(vapply(groups, function(ii) mean(log(m$library_size[ii])), numeric(1))))
  exposure <- m$library_size / reference
  rate <- colMeans(do.call(rbind, lapply(groups, function(ii) colMeans(X[ii,,drop=FALSE] / exposure[ii]))))
  mu <- exposure %o% pmax(rate, 1e-8)
  raw <- colMeans(do.call(rbind, lapply(groups, function(ii) {
    pmax(colMeans((X[ii,,drop=FALSE] - mu[ii,,drop=FALSE])^2 - mu[ii,,drop=FALSE]) /
           pmax(colMeans(mu[ii,,drop=FALSE]^2), 1e-8), 0)
  })))
  target <- median(raw[is.finite(raw) & raw > 0])
  if (!is.finite(target)) target <- 0.1
  phi <- pmin(pmax(0.5 * raw + 0.5 * target, 1e-8), 100)
  Y <- (X - mu) / sqrt(pmax(mu + sweep(mu^2, 2, phi, '*'), 1e-8))
  bound <- sqrt(nrow(X) / 30)
  pmin(pmax(Y, -bound), bound)
}

fit_regift_input <- function(x, dataset = c('kang_ifnb', 'parse_tnf', 'hpap_t1d'),
                             max_iter = 100L, threads = 8L) {
  dataset <- match.arg(dataset)
  stopifnot(as.character(utils::packageVersion('ReGIFT')) == '0.1.0',
            nrow(x$counts) == nrow(x$meta), ncol(x$counts) == length(x$genes))
  Y <- regift_input_working(x, dataset)
  set.seed(if (dataset == 'parse_tnf') 20260905L else 20260909L)
  args <- list(Y = Y, meta = x$meta, contrasts = x$contrasts,
               K = if (dataset == 'hpap_t1d') 10L else 5L, H = 2L,
               lambda_Delta = 3, kappa = 1.345, tol = 1e-6,
               z_backend = 'blocked', gene_block = 256L,
               threads = as.integer(threads), svd_backend = 'randomized')
  probe <- do.call(ReGIFT::regift_fit, c(args, list(lambda_B = 0,
                    max_iter = 1L, state_anchor_shrinkage = FALSE)))
  do.call(ReGIFT::regift_fit, c(args, list(lambda_B = probe$lambda_B_max / 32,
                    max_iter = as.integer(max_iter), state_anchor_shrinkage = TRUE)))
}
