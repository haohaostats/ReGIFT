#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Dense correctness-preserving kernel for the cell-state update. Arrays use
// R column-major indexing: B[p,K,Q] and Delta[S,p,K,Q].
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_update_z(
    const Rcpp::NumericMatrix& Y_r,
    const Rcpp::NumericMatrix& nuisance_r,
    const Rcpp::NumericMatrix& L0_r,
    const Rcpp::NumericVector& B,
    const Rcpp::NumericVector& Delta,
    const Rcpp::NumericMatrix& D_r,
    const Rcpp::IntegerVector& donor,
    const double ridge = 1e-6,
    const int threads = 1) {

  const int n = Y_r.nrow();
  const int p = Y_r.ncol();
  const int K = L0_r.ncol();
  const int Q = D_r.ncol();
  const Rcpp::IntegerVector bdim = B.attr("dim");
  const Rcpp::IntegerVector ddim = Delta.attr("dim");
  if (bdim.size() != 3 || bdim[0] != p || bdim[1] != K || bdim[2] != Q)
    Rcpp::stop("B must have dimensions p x K x Q.");
  if (ddim.size() != 4 || ddim[1] != p || ddim[2] != K || ddim[3] != Q)
    Rcpp::stop("Delta must have dimensions S x p x K x Q.");
  const int S = ddim[0];
  if (nuisance_r.nrow() != n || nuisance_r.ncol() != p ||
      D_r.nrow() != n || donor.size() != n)
    Rcpp::stop("Cell-level inputs have incompatible dimensions.");

  arma::mat Y(const_cast<double*>(Y_r.begin()), n, p, false, true);
  arma::mat nuisance(const_cast<double*>(nuisance_r.begin()), n, p, false, true);
  arma::mat L0(const_cast<double*>(L0_r.begin()), p, K, false, true);
  arma::mat D(const_cast<double*>(D_r.begin()), n, Q, false, true);
  Rcpp::NumericMatrix out(n, K);

#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#pragma omp parallel for schedule(static) if(threads > 1)
#endif
  for (int i = 0; i < n; ++i) {
    const int s = donor[i] - 1;
    if (s < 0 || s >= S) continue;
    arma::mat loading = L0;
    for (int q = 0; q < Q; ++q) {
      const double dq = D(i, q);
      if (dq == 0.0) continue;
      for (int k = 0; k < K; ++k) {
        for (int g = 0; g < p; ++g) {
          const R_xlen_t ib = g + static_cast<R_xlen_t>(p) * (k + K * q);
          const R_xlen_t id = s + static_cast<R_xlen_t>(S) *
            (g + static_cast<R_xlen_t>(p) * (k + K * q));
          loading(g, k) += dq * (B[ib] + Delta[id]);
        }
      }
    }
    arma::mat gram = loading.t() * loading;
    gram.diag() += ridge;
    arma::vec rhs = loading.t() * (Y.row(i) - nuisance.row(i)).t();
    arma::vec zi;
    bool ok = arma::solve(zi, gram, rhs, arma::solve_opts::likely_sympd);
    if (!ok) zi = arma::pinv(gram) * rhs;
    for (int k = 0; k < K; ++k) out(i, k) = zi[k];
  }
  return out;
}

// Gene-blocked Z update. Unlike cpp_update_z, this never allocates a p x K
// loading matrix for a cell. It accumulates the K x K Gram matrix and K-vector
// right-hand side one gene block at a time and is algebraically identical.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_update_z_blocked(
    const Rcpp::NumericMatrix& Y,
    const Rcpp::NumericMatrix& nuisance,
    const Rcpp::NumericMatrix& L0,
    const Rcpp::NumericVector& B,
    const Rcpp::NumericVector& Delta,
    const Rcpp::NumericMatrix& D,
    const Rcpp::IntegerVector& donor,
    const int gene_block = 256,
    const double ridge = 1e-6,
    const int threads = 1) {
  const int n=Y.nrow(), p=Y.ncol(), K=L0.ncol(), Q=D.ncol();
  const Rcpp::IntegerVector bd=B.attr("dim"), dd=Delta.attr("dim");
  if (nuisance.nrow()!=n || nuisance.ncol()!=p || L0.nrow()!=p ||
      D.nrow()!=n || donor.size()!=n || bd.size()!=3 || dd.size()!=4)
    Rcpp::stop("Blocked Z inputs have incompatible dimensions.");
  const int S=dd[0], block=std::max(1,gene_block);
  Rcpp::NumericMatrix out(n,K);
#ifdef _OPENMP
  if (threads>0) omp_set_num_threads(threads);
#pragma omp parallel for schedule(static) if(threads > 1)
#endif
  for (int i=0;i<n;++i) {
    const int s=donor[i]-1;
    if (s<0 || s>=S) continue;
    arma::mat gram(K,K,arma::fill::zeros);
    arma::vec rhs(K,arma::fill::zeros), load(K);
    for (int first=0;first<p;first+=block) {
      const int last=std::min(p,first+block);
      for (int g=first;g<last;++g) {
        for (int k=0;k<K;++k) {
          double value=L0(g,k);
          for (int q=0;q<Q;++q) {
            const R_xlen_t ib=g+static_cast<R_xlen_t>(p)*(k+K*q);
            const R_xlen_t id=s+static_cast<R_xlen_t>(S)*(g+static_cast<R_xlen_t>(p)*(k+K*q));
            value += D(i,q)*(B[ib]+Delta[id]);
          }
          load[k]=value;
        }
        gram += load*load.t();
        rhs += load*(Y(i,g)-nuisance(i,g));
      }
    }
    gram.diag() += ridge;
    arma::vec zi;
    bool ok=arma::solve(zi,gram,rhs,arma::solve_opts::likely_sympd);
    if(!ok) zi=arma::pinv(gram)*rhs;
    for(int k=0;k<K;++k) out(i,k)=zi[k];
  }
  return out;
}
