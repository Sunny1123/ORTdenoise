// // [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// // [[Rcpp::plugins(cpp11)]]
//
// #include <RcppArmadillo.h>
// #include <RcppParallel.h>
//
// using namespace Rcpp;
// using namespace RcppParallel;
//
// // Compute integral image
// arma::mat integral_image(const arma::mat &img) {
//   arma::mat S = arma::zeros<arma::mat>(img.n_rows, img.n_cols);
//   S.row(0) = arma::cumsum(img.row(0), 1);
//   for (size_t i=1; i<img.n_rows; i++) {
//     S.row(i) = S.row(i-1) + arma::cumsum(img.row(i),1);
//   }
//   return S;
// }
//
// // Sum over rectangle using integral image
// inline double rect_sum(const arma::mat &S, int x1, int y1, int x2, int y2) {
//   double total = S(x2,y2);
//   if (x1>0) total -= S(x1-1,y2);
//   if (y1>0) total -= S(x2,y1-1);
//   if (x1>0 && y1>0) total += S(x1-1,y1-1);
//   return total;
// }
//
// // Worker for parallel SSIM over channels
// struct SSIMWorker : public Worker {
//   const arma::cube &img1;
//   const arma::cube &img2;
//   int win_size;
//   double C1, C2;
//   arma::vec &ssim_channel;
//
//   SSIMWorker(const arma::cube &img1_,
//              const arma::cube &img2_,
//              int win_size_,
//              double K1, double K2,
//              arma::vec &ssim_channel_)
//     : img1(img1_), img2(img2_), win_size(win_size_),
//       C1(K1*K1), C2(K2*K2), ssim_channel(ssim_channel_) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     int H = img1.n_rows;
//     int W = img1.n_cols;
//     int pad = win_size/2;
//     int area = win_size * win_size;
//
//     for (size_t c=begin; c<end; c++) {
//       arma::mat x = img1.slice(c);
//       arma::mat y = img2.slice(c);
//
//       // Compute integral images
//       arma::mat Sx = integral_image(x);
//       arma::mat Sy = integral_image(y);
//       arma::mat Sx2 = integral_image(arma::square(x));
//       arma::mat Sy2 = integral_image(arma::square(y));
//       arma::mat Sxy = integral_image(x % y);
//
//       arma::mat ssim_map = arma::zeros<arma::mat>(H,W);
//
//       // Compute SSIM map
//       for (int i=pad; i<H-pad; i++) {
//         for (int j=pad; j<W-pad; j++) {
//           int x1=i-pad, y1=j-pad, x2=i+pad, y2=j+pad;
//
//           double sum_x  = rect_sum(Sx, x1,y1,x2,y2);
//           double sum_y  = rect_sum(Sy, x1,y1,x2,y2);
//           double sum_x2 = rect_sum(Sx2,x1,y1,x2,y2);
//           double sum_y2 = rect_sum(Sy2,x1,y1,x2,y2);
//           double sum_xy = rect_sum(Sxy,x1,y1,x2,y2);
//
//           double mu_x = sum_x / area;
//           double mu_y = sum_y / area;
//           double sigma_x2 = sum_x2/area - mu_x*mu_x;
//           double sigma_y2 = sum_y2/area - mu_y*mu_y;
//           double sigma_xy = sum_xy/area - mu_x*mu_y;
//
//           ssim_map(i,j) = ((2*mu_x*mu_y + C1)*(2*sigma_xy + C2)) /
//             ((mu_x*mu_x + mu_y*mu_y + C1)*(sigma_x2 + sigma_y2 + C2));
//         }
//       }
//
//       // average SSIM ignoring border padding
//       ssim_channel(c) = arma::mean(arma::vectorise(ssim_map.rows(pad,H-pad-1).cols(pad,W-pad-1)));
//     }
//   }
// };
//
// // [[Rcpp::export]]
// double ssim(arma::cube img1, arma::cube img2,
//             double K1=0.01, double K2=0.03, int win_size=7) {
//
//   if (img1.n_rows != img2.n_rows || img1.n_cols != img2.n_cols || img1.n_slices != img2.n_slices)
//     stop("Dimensions must match");
//
//   int C = img1.n_slices;
//   arma::vec ssim_channel(C, arma::fill::zeros);
//
//   SSIMWorker worker(img1, img2, win_size, K1, K2, ssim_channel);
//   parallelFor(0,C,worker);
//
//   return arma::mean(ssim_channel);
// }


// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// Compute integral image
arma::mat integral_image(const arma::mat &img) {
  arma::mat S = arma::zeros<arma::mat>(img.n_rows, img.n_cols);
  S.row(0) = arma::cumsum(img.row(0), 1);
  for (size_t i=1; i<img.n_rows; i++) {
    S.row(i) = S.row(i-1) + arma::cumsum(img.row(i),1);
  }
  return S;
}

// Sum over rectangle using integral image
inline double rect_sum(const arma::mat &S, int x1, int y1, int x2, int y2) {
  double total = S(x2,y2);
  if (x1>0) total -= S(x1-1,y2);
  if (y1>0) total -= S(x2,y1-1);
  if (x1>0 && y1>0) total += S(x1-1,y1-1);
  return total;
}

// Worker for parallel SSIM over channels
struct SSIMWorker : public Worker {
  const arma::cube &img1;
  const arma::cube &img2;
  int win_size;
  double C1, C2;
  arma::vec &ssim_channel;

  SSIMWorker(const arma::cube &img1_,
             const arma::cube &img2_,
             int win_size_,
             double K1, double K2,
             arma::vec &ssim_channel_)
    : img1(img1_), img2(img2_), win_size(win_size_),
      C1(K1*K1), C2(K2*K2), ssim_channel(ssim_channel_) {}

  void operator()(std::size_t begin, std::size_t end) {
    int H = img1.n_rows;
    int W = img1.n_cols;
    int pad = win_size/2;
    int area = win_size * win_size;

    for (size_t c=begin; c<end; c++) {
      arma::mat x = img1.slice(c);
      arma::mat y = img2.slice(c);

      arma::mat Sx = integral_image(x);
      arma::mat Sy = integral_image(y);
      arma::mat Sx2 = integral_image(arma::square(x));
      arma::mat Sy2 = integral_image(arma::square(y));
      arma::mat Sxy = integral_image(x % y);

      arma::mat ssim_map = arma::zeros<arma::mat>(H,W);

      for (int i=pad; i<H-pad; i++) {
        for (int j=pad; j<W-pad; j++) {
          int x1=i-pad, y1=j-pad, x2=i+pad, y2=j+pad;

          double sum_x  = rect_sum(Sx, x1,y1,x2,y2);
          double sum_y  = rect_sum(Sy, x1,y1,x2,y2);
          double sum_x2 = rect_sum(Sx2,x1,y1,x2,y2);
          double sum_y2 = rect_sum(Sy2,x1,y1,x2,y2);
          double sum_xy = rect_sum(Sxy,x1,y1,x2,y2);

          double mu_x = sum_x / area;
          double mu_y = sum_y / area;
          double sigma_x2 = sum_x2/area - mu_x*mu_x;
          double sigma_y2 = sum_y2/area - mu_y*mu_y;
          double sigma_xy = sum_xy/area - mu_x*mu_y;

          ssim_map(i,j) = ((2*mu_x*mu_y + C1)*(2*sigma_xy + C2)) /
            ((mu_x*mu_x + mu_y*mu_y + C1)*(sigma_x2 + sigma_y2 + C2));
        }
      }
      ssim_channel(c) = arma::mean(arma::vectorise(
        ssim_map.rows(pad,H-pad-1).cols(pad,W-pad-1)
      ));
    }
  }
};

// [[Rcpp::export]]
double ssim(SEXP img1_, SEXP img2_,
            double K1=0.01, double K2=0.03, int win_size=7) {

  arma::cube img1, img2;

  // handle grayscale (matrix) → convert to cube with 1 slice
  if (Rf_isMatrix(img1_) && Rf_isMatrix(img2_)) {
    arma::mat m1 = as<arma::mat>(img1_);
    arma::mat m2 = as<arma::mat>(img2_);
    if (m1.n_rows != m2.n_rows || m1.n_cols != m2.n_cols)
      stop("Dimensions must match");
    img1.set_size(m1.n_rows, m1.n_cols, 1);
    img2.set_size(m2.n_rows, m2.n_cols, 1);
    img1.slice(0) = m1;
    img2.slice(0) = m2;
  } else {
    img1 = as<arma::cube>(img1_);
    img2 = as<arma::cube>(img2_);
    if (img1.n_rows != img2.n_rows ||
        img1.n_cols != img2.n_cols ||
        img1.n_slices != img2.n_slices)
      stop("Dimensions must match");
  }

  int C = img1.n_slices;
  arma::vec ssim_channel(C, arma::fill::zeros);

  SSIMWorker worker(img1, img2, win_size, K1, K2, ssim_channel);
  parallelFor(0,C,worker);

  return arma::mean(ssim_channel);
}
