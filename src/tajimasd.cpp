#include <Rcpp.h>

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

// https://gallery.rcpp.org/articles/parallel-distance-matrix/

template <typename InputIterator1, typename InputIterator2>
inline double tajimas_d(InputIterator1 begin1, InputIterator1 end1,
                        InputIterator2 begin2) {
  // we are trying to calculate "d" like here:
  // https://en.wikipedia.org/wiki/Tajima%27s_D

  // the goal is to find the number of bases that differ between the two seqs
  // and then divide by the effective number of comparisons
  // so the denominator should be 1 when there is no missing data
  // and 0.8 when there is 20% incomparably positons, etc

  // the code is very similar to the identity distance calculation
  int numer = 0;
  int denom = 0;
  int denom_max = 0;

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;

  // loop over every combination of items
  while (it1 != end1) {
    // dereference the iterators and increment them
    int d1 = *it1++;
    int d2 = *it2++;

    denom_max += 1;

    if (d1 != -9 && d2 != -9) {
      denom += 1;
      if (d1 != d2) {
        numer += 1;
      }
    }
  }

  double k;
  if (denom == 0) {
    // guard against div by zero
    // which means there were no comparable positions
    k = 0.0;
  } else {
    // have to static_cast<double> in order to get floating point division out of two ints
    k = static_cast<double>(numer)/(static_cast<double>(denom)/static_cast<double>(denom_max));
  }

  return k;
}

//'
//' rcpp_tajimas_d
//'
//' @name rcpp_tajimas_d
//' @param mat The input matrix (k markers times n individuals).
//' @return rmat The output matrix (n individuals times n)
//'
//' @export
// [[Rcpp::export]]
double rcpp_tajimas_d(Rcpp::IntegerMatrix mat) {
  // allocate the return matrix
  Rcpp::NumericMatrix rmat(mat.ncol(), mat.ncol());

  for (int i = 0; i < rmat.nrow(); i++) {
    for (int j = 0; j < i; j++) {
      // get the iterable columns we operate on
      Rcpp::IntegerMatrix::Column col1 = mat.column(i);
      Rcpp::IntegerMatrix::Column col2 = mat.column(j);
      rmat(j, i) = tajimas_d(col1.begin(), col1.end(), col2.begin());
      rmat(i, j) = rmat(j, i);
    }
  }

  rownames(rmat) = colnames(mat);
  colnames(rmat) = colnames(mat);

  // calculate khat, which is like pi
  double khat = static_cast<double>(Rcpp::sum(rmat))/(mat.ncol() * (mat.ncol() - 1));

  // calculate the number of segregating sites, `S`
  // and divide by a1= \sum_{i=1}^{n-1} 1/i
  int S = 0;

  // find each row where at least one individual has a different allele
  for (int i = 0; i < mat.nrow(); i++) {
    for (int j = 0; j < mat.ncol(); j++) {
      if (mat(i, j) != 0 && mat(i, j) != -9) {
        S += 1;
        break;
      }
    }
  }

  // calculate denominator a1
  // and axillary a2
  double a1 = 0;
  double a2 = 0;
  for (int i = 1; i < mat.ncol(); i++) {
    a1 += 1/static_cast<double>(i);
    a2 += 1/std::pow(static_cast<double>(i), 2);
  }

  double M = S/a1;

  double d = khat - M;

  double n = static_cast<double>(mat.ncol());
  double b1 = (n + 1.0)/(3.0 * (n - 1.0));
  double b2 = 2.0 * (std::pow(n, 2) + n + 3.0)/(9.0 * n * (n - 1.0));
  double c1 = b1 - 1/a1;
  double c2 = b2 - (n + 2.0)/(a1 * n) + (a2/std::pow(a1, 2));
  double e1 = c1/a1;
  double e2 = c2/(std::pow(a1, 2) + a2);

  double denom = std::sqrt(e1 * static_cast<double>(S) +
                           e2 * static_cast<double>(S) * (static_cast<double>(S - 1.0)));

  double D = d/denom;
  return D;
}

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

struct TajimasD : public RcppParallel::Worker {
  // input mat
  const RcppParallel::RMatrix<int> mat;

  // output mat
  RcppParallel::RMatrix<double> rmat;

  // some weird boilerplate to make it work?
  TajimasD(const Rcpp::IntegerMatrix mat, Rcpp::NumericMatrix rmat)
    : mat(mat), rmat(rmat) {}

  // the fun call operator will only operate on _certain columns_
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < i; j ++) {
        RcppParallel::RMatrix<int>::Column col1 = mat.column(i);
        RcppParallel::RMatrix<int>::Column col2 = mat.column(j);
        rmat(j, i) = tajimas_d(col1.begin(), col1.end(), col2.begin());
        rmat(i, j) = rmat(j, i);
      }
    }
  }
};

//'
//' rcpp_parallel_tajimas_d
//'
//' @name rcpp_parallel_tajimas_d
//' @param mat The input matrix (k markers times n individuals).
//' @return rmat The output vector
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_parallel_tajimas_d(Rcpp::IntegerMatrix mat) {
  Rcpp::NumericMatrix rmat(mat.ncol(), mat.ncol());

  TajimasD tajimasD(mat, rmat);

  RcppParallel::parallelFor(0, mat.ncol(), tajimasD);

  rownames(rmat) = colnames(mat);
  colnames(rmat) = colnames(mat);

  // calculate khat, which is like pi
  double khat = static_cast<double>(Rcpp::sum(rmat))/(mat.ncol() * (mat.ncol() - 1));

  // calculate the number of segregating sites, `S`
  // and divide by a1= \sum_{i=1}^{n-1} 1/i
  int S = 0;

  // find each row where at least one individual has a different allele
  for (int i = 0; i < mat.nrow(); i++) {
    for (int j = 0; j < mat.ncol(); j++) {
      if (mat(i, j) != 0 && mat(i, j) != -9) {
        S += 1;
        break;
      }
    }
  }

  // calculate denominator a1
  // and axillary a2
  double a1 = 0;
  double a2 = 0;
  for (int i = 1; i < mat.ncol(); i++) {
    a1 += 1/static_cast<double>(i);
    a2 += 1/std::pow(static_cast<double>(i), 2);
  }

  double M = S/a1;

  double d = khat - M;

  double n = static_cast<double>(mat.ncol());
  double b1 = (n + 1.0)/(3.0 * (n - 1.0));
  double b2 = 2.0 * (std::pow(n, 2) + n + 3.0)/(9.0 * n * (n - 1.0));
  double c1 = b1 - 1/a1;
  double c2 = b2 - (n + 2.0)/(a1 * n) + (a2/std::pow(a1, 2));
  double e1 = c1/a1;
  double e2 = c2/(std::pow(a1, 2) + a2);

  double denom = std::sqrt(e1 * static_cast<double>(S) +
                           e2 * static_cast<double>(S) * (static_cast<double>(S - 1.0)));

  double D;

  if (denom == 0.0) {
    D = 0.0;
  } else {
    D = d/denom;
  }

  double pi = khat/mat.nrow() * (mat.ncol()/(mat.ncol() - 1.0));

  Rcpp::NumericVector outvec = Rcpp::NumericVector::create(
    Rcpp::Named("pi", pi),
    Rcpp::Named("khat", khat),
    Rcpp::Named("S", S),
    Rcpp::Named("N", mat.nrow()),
    Rcpp::Named("M", M),
    Rcpp::Named("d", d),
    Rcpp::Named("D", D));
  return outvec;
}
