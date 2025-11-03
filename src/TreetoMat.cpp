#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

using namespace Rcpp;
using namespace RcppParallel;

// ==========================================================
// Helpers
// ==========================================================
inline int ravel_index(const std::vector<int>& coords,
                       const std::vector<int>& dims) {
  int idx = 0, mult = 1;
  for (size_t k = 0; k < dims.size(); ++k) {
    idx += coords[k] * mult;
    mult *= dims[k];
  }
  return idx;
}

// Generate all combinations of offsets in [-r, r] for N dimensions
void generate_offsets(int N, int r, std::vector<std::vector<int>>& offsets) {
  std::vector<int> current(N, -r);
  bool done = false;
  while (!done) {
    offsets.push_back(current);
    for (int i = 0; i < N; ++i) {
      current[i]++;
      if (current[i] <= r) break;
      current[i] = -r;
      if (i == N-1) done = true;
    }
  }
}

// ==========================================================
// Lookup: precomputed exp table for speed
// ==========================================================
// struct Lookup {
//   std::vector<double> exp_table;
//   double step_ratio;
//   int n_ratio;
//   int smband;
//
//   Lookup(double maxRatio, int nSteps, int smband_)
//     : step_ratio(nSteps > 1 ? maxRatio / (nSteps - 1) : 1.0),
//       n_ratio(std::max(1, nSteps)),
//       smband(smband_) {
//     exp_table.resize(n_ratio);
//     for (int i = 0; i < n_ratio; ++i)
//       exp_table[i] = std::exp(-(i * step_ratio));
//   }
//
//   inline double getExp(double ratio) const {
//     if (ratio <= 0.0) return 1.0;
//     if (step_ratio <= 0.0) return std::exp(-ratio);
//     int idx = (int)(ratio / step_ratio);
//     if (idx < 0) idx = 0;
//     if (idx >= n_ratio) idx = n_ratio - 1;
//     return exp_table[idx];
//   }
// };
struct Lookup {
  std::vector<double> exp_table;
  double step_ratio;
  int n_ratio;
  int smband;
  double factor;   // new: multiplicative factor for exponent

  Lookup(double maxRatio, int nSteps, int smband_, double factor_ = 100.0)
    : step_ratio(nSteps > 1 ? maxRatio / (nSteps - 1) : 1.0),
      n_ratio(std::max(1, nSteps)),
      smband(smband_),
      factor(factor_)
  {
    exp_table.resize(n_ratio);
    for (int i = 0; i < n_ratio; ++i)
      exp_table[i] = std::exp(-factor * (i * step_ratio));
  }

  inline double getExp(double ratio) const {
    if (ratio <= 0.0) return 1.0;
    if (step_ratio <= 0.0) return std::exp(-factor * ratio);
    int idx = (int)(ratio / step_ratio);
    if (idx < 0) idx = 0;
    if (idx >= n_ratio) idx = n_ratio - 1;
    return exp_table[idx];
  }
};
// ==========================================================
// PatchDatabase: N-dimensional patch storage
// ==========================================================
struct PatchDatabase {
  std::vector<double> patches;
  std::vector<double> norms;
  int patchRadius;
  std::vector<int> dims;
  int patchSize;
  int N;

  PatchDatabase(const RVector<double>& img,
                const std::vector<int>& dims_vec,
                int patchRadius_)
    : patchRadius(patchRadius_), dims(dims_vec) {

    N = 1;
    for (size_t i = 0; i < dims.size(); ++i) N *= dims[i];

    patchSize = 1;
    for (size_t i = 0; i < dims.size(); ++i)
      patchSize *= (2 * patchRadius + 1);

    patches.assign((size_t)N * patchSize, 0.0);
    norms.assign(N, 0.0);

    std::vector<std::vector<int>> offsets;
    generate_offsets(dims.size(), patchRadius, offsets);

    std::vector<int> coords(dims.size());
    for (int lin = 0; lin < N; ++lin) {
      int rem = lin;
      for (size_t d = 0; d < dims.size(); ++d) {
        coords[d] = rem % dims[d];
        rem /= dims[d];
      }

      int k = 0;
      double sumSq = 0.0;
      for (const auto& off : offsets) {
        std::vector<int> pn = coords;
        for (size_t d = 0; d < dims.size(); ++d) {
          pn[d] += off[d];
          if (pn[d] < 0) pn[d] = 0;
          if (pn[d] >= dims[d]) pn[d] = dims[d] - 1;
        }
        int lin_p = ravel_index(pn, dims);
        double v = img[lin_p];
        patches[(size_t)lin * patchSize + k] = v;
        sumSq += v * v;
        ++k;
      }
      norms[lin] = sumSq;
    }
  }

  inline double dot_patch(int i, int j) const {
    const double* Pi = &patches[(size_t)i * patchSize];
    const double* Pj = &patches[(size_t)j * patchSize];
    double s = 0.0;
    for (int k = 0; k < patchSize; ++k) s += Pi[k] * Pj[k];
    return s;
  }

  inline double distance2(int i, int j) const {
    double d2 = norms[i] + norms[j] - 2.0 * dot_patch(i, j);
    return (d2 < 0.0 ? 0.0 : d2);
  }
};

// ==========================================================
// PixelWorker: N-dimensional neighbor processing
// ==========================================================
struct PixelWorker : public RcppParallel::Worker {
  const RVector<double> img;
  const std::vector<std::vector<int>>& pixel_coords;
  const std::vector<int>& pixel_nodes;
  const std::vector<int>& dims;
  const Lookup& lookup;
  const PatchDatabase& pdb;
  RVector<double> output;
  double h2;

  PixelWorker(const RVector<double>& img_,
              const std::vector<std::vector<int>>& pixel_coords_,
              const std::vector<int>& pixel_nodes_,
              const std::vector<int>& dims_,
              const Lookup& lookup_,
              const PatchDatabase& pdb_,
              NumericVector& output_,
              double h2_)
    : img(img_), pixel_coords(pixel_coords_), pixel_nodes(pixel_nodes_),
      dims(dims_), lookup(lookup_), pdb(pdb_), output(output_), h2(h2_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const double eps = 1e-14;
    int Ndim = dims.size();
    std::vector<std::vector<int>> offsets;
    generate_offsets(Ndim, lookup.smband, offsets);

    for (std::size_t lin = begin; lin < end; ++lin) {
      const auto& ci = pixel_coords[lin];
      int node_id = pixel_nodes[lin];
      double centerVal = img[lin];

      double temp = 0.0;
      double count = 0.0;

      for (const auto& off : offsets) {
        std::vector<int> neigh = ci;
        for (int d = 0; d < Ndim; ++d) {
          neigh[d] += off[d];
          if (neigh[d] < 0) neigh[d] = 0;
          if (neigh[d] >= dims[d]) neigh[d] = dims[d] - 1;
        }

        int lin_j = ravel_index(neigh, dims);
        if ((int)lin_j == (int)lin) continue;
        if (pixel_nodes[lin_j] != node_id) continue;

        double pd2 = pdb.distance2((int)lin, lin_j);
        if (pd2 <= eps) continue;

        double ratio = pd2 / h2;
        double wt = lookup.getExp(ratio);
        if (wt < 1e-12) continue;

        temp += img[lin_j] * wt;
        count += wt;
      }

      output[lin] = (count > 0.0 ? temp / count : centerVal);
    }
  }
};

// ==========================================================
// Exported function: N-dimensional support
// ==========================================================
// [[Rcpp::export]]
NumericVector TreetomatPixelSIMD(List tree,
                                 NumericVector img,
                                 IntegerVector dims,
                                 double smband,
                                 int nSteps = 10000) {
  std::vector<int> dims_vec = as<std::vector<int>>(dims);
  int N = std::accumulate(dims_vec.begin(), dims_vec.end(), 1, std::multiplies<int>());
  if (img.size() != N) stop("Image size does not match dims");

  NumericVector output(N);
  const int patchRadius = 1;

  double imgMin = *std::min_element(img.begin(), img.end());
  double imgMax = *std::max_element(img.begin(), img.end());

  double h = 0.6 * (imgMax - imgMin) * std::sqrt((double)std::pow((2*patchRadius+1), dims_vec.size()));
  if (h <= 0.0) h = 1.0;
  double h2 = h*h;
  int scalefactor=5;
  double maxPatchDist2 = std::pow((2*patchRadius+1), dims_vec.size()) * (imgMax - imgMin) * (imgMax - imgMin);
  double maxRatio = maxPatchDist2 / scalefactor;
  if (maxRatio <= 0.0) maxRatio = 1.0;

  Lookup lookup(maxRatio, nSteps, (int)std::ceil(smband),scalefactor);

  std::vector<std::vector<int>> pixel_coords(N);
  std::vector<int> pixel_nodes(N, -1);

  for (int nodeIdx = 0; nodeIdx < tree.size(); ++nodeIdx) {
    List node = tree[nodeIdx];
    IntegerMatrix node_coords = node["coord"];
    for (int i = 0; i < node_coords.nrow(); ++i) {
      std::vector<int> coords(dims_vec.size());
      for (size_t d = 0; d < dims_vec.size(); ++d)
        coords[d] = node_coords(i, d+1) - 1;
      int lin_idx = ravel_index(coords, dims_vec);
      pixel_coords[lin_idx] = coords;
      pixel_nodes[lin_idx] = nodeIdx;
    }
  }

  // fallback for unassigned pixels
  for (int idx = 0; idx < N; ++idx) {
    if (pixel_nodes[idx] == -1) {
      if (pixel_coords[idx].empty()) {
        std::vector<int> coords(dims_vec.size());
        int rem = idx;
        for (size_t d = 0; d < dims_vec.size(); ++d) {
          coords[d] = rem % dims_vec[d];
          rem /= dims_vec[d];
        }
        pixel_coords[idx] = coords;
      }
      pixel_nodes[idx] = idx;
    }
  }

  RVector<double> imgR(img);
  PatchDatabase pdb(imgR, dims_vec, patchRadius);

  PixelWorker worker(imgR, pixel_coords, pixel_nodes, dims_vec, lookup, pdb, output, h2);
  parallelFor(0, N, worker);

  return output;
}
