#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <numeric>

using namespace Rcpp;
using namespace RcppParallel;

// ---------------- Node structure ----------------
struct Node {
  arma::vec intensity;
  arma::mat coord;
  bool leaf;
};

// ---------------- Loss function ----------------
double loss_eval(const Node &node, const arma::vec &w) {
  const arma::vec &y = node.intensity;
  const arma::mat &X = node.coord;
  int n = y.n_elem;
  int p = w.n_elem;

  if (n <= 1) return 0.0;

  arma::vec w_norm = w / arma::norm(w);
  arma::vec xmod = X * w_norm;
  arma::uvec idx_left = arma::find(xmod > 0);
  arma::uvec idx_right = arma::find(xmod <= 0);

  int n1 = idx_left.n_elem;
  int n2 = idx_right.n_elem;

  if (n1 <= 1 || n2 <= 1) return 0.0;

  double mean1 = arma::mean(y.elem(idx_left));
  double mean2 = arma::mean(y.elem(idx_right));
  double r1 = double(n1) / n;
  double r2 = double(n2) / n;

  return r1 * r2 * std::pow(mean1 - mean2, 2);
}

// ---------------- Nelder-Mead optimizer ----------------
arma::vec optimize_w_max(const Node &node, double tol=1e-5, int maxit=10000) {
  int p = node.coord.n_cols;

  // Initial guess
  arma::vec w_init(p, arma::fill::ones);
  if (p > 1) {
    arma::mat sub = node.coord.cols(1, p - 1);
    double col_mean = arma::accu(sub) / sub.n_elem;
    w_init[0] = - (p - 1) * col_mean;
    for (size_t i = 1; i < p; i++) w_init[i] = 1.0;
  }
  w_init = w_init / arma::norm(w_init);

  // Simplex initialization
  std::vector<arma::vec> simplex(p+1, w_init);
  double scale = 0.05;
  for (int i = 1; i <= p; i++) simplex[i][i-1] += scale;

  // Track global best
  arma::vec best_w = w_init;
  double best_loss = -loss_eval(node, w_init); // store negative since we are minimizing
  int nsteps = 0;

  while (nsteps < maxit) {
    std::vector<double> fvals(p+1);
    for (int i = 0; i <= p; i++) fvals[i] = -loss_eval(node, simplex[i]); // negative

    // Update global best
    for (int i = 0; i <= p; i++) {
      if (fvals[i] < best_loss) { // remember fvals is negative
        best_loss = fvals[i];
        best_w = simplex[i];
      }
    }

    // Order simplex
    std::vector<int> idx(p+1);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return fvals[a] < fvals[b]; });

    arma::vec x_best = simplex[idx[0]];
    arma::vec x_worst = simplex[idx[p]];
    arma::vec x_second_worst = simplex[idx[p-1]];

    // Centroid
    arma::vec centroid = arma::zeros(p);
    for (int i = 0; i < p; i++) centroid += simplex[idx[i]];
    centroid /= p;

    double alpha = 1.0;
    arma::vec xr = centroid + alpha * (centroid - x_worst);
    double fr = -loss_eval(node, xr);

    if (fr < fvals[idx[0]]) {
      double gamma = 2.0;
      arma::vec xe = centroid + gamma * (xr - centroid);
      double fe = -loss_eval(node, xe);
      if (fe < fr) simplex[idx[p]] = xe;
      else simplex[idx[p]] = xr;
    } else if (fr < fvals[idx[p-1]]) {
      simplex[idx[p]] = xr;
    } else {
      double rho = 0.5;
      arma::vec xc;
      if (fr < fvals[idx[p]]) xc = centroid + rho * (xr - centroid);
      else xc = centroid + rho * (x_worst - centroid);
      double fc = -loss_eval(node, xc);
      if (fc < fvals[idx[p]]) simplex[idx[p]] = xc;
      else {
        for (int i = 1; i <= p; i++)
          simplex[idx[i]] = x_best + 0.5*(simplex[idx[i]] - x_best);
      }
    }

    // Ensure valid split
    arma::vec xmod = node.coord * (best_w / arma::norm(best_w));
    arma::uvec idx1 = arma::find(xmod > 0);
    arma::uvec idx2 = arma::find(xmod <= 0);
    if (idx1.n_elem == 0 || idx2.n_elem == 0) {
      for (int j = 1; j < p; j++) best_w[j] = R::runif(-1.0, 1.0);
      best_w[0] = -arma::mean(node.coord * best_w);
    }

    // Check convergence
    double fstd = 0;
    for (double v : fvals) fstd += std::pow(v - fvals[idx[0]], 2);
    fstd = std::sqrt(fstd / (p+1));
    if (fstd < tol) break;

    nsteps++;
  }

  return best_w / arma::norm(best_w);
}

// ---------------- Get children ----------------
std::vector<Node> getchild(const Node &node, int minsize, double minlossimprove) {
  if (node.leaf) return {node};

  arma::vec w_opti = optimize_w_max(node,minlossimprove*0.1);
  arma::vec xmod = node.coord * w_opti;
  arma::uvec idx1 = arma::find(xmod > 0);
  arma::uvec idx2 = arma::find(xmod <= 0);

  Node child_left = {node.intensity.elem(idx1), node.coord.rows(idx1), false};
  Node child_right = {node.intensity.elem(idx2), node.coord.rows(idx2), false};

  double lsimp = 0.0;
  if (idx1.n_elem > 1 && idx2.n_elem > 1) lsimp = loss_eval(node, w_opti);

  if (lsimp > minlossimprove &&
      child_left.intensity.n_elem > minsize &&
      child_right.intensity.n_elem > minsize) {
    return {child_left, child_right};
  } else {
    Node leaf_node = node;
    leaf_node.leaf = true;
    return {leaf_node};
  }
}

// ---------------- Parallel worker ----------------
struct NodeWorker : public Worker {
  std::vector<Node>& nodes_in;
  std::vector<std::vector<Node>>& nodes_out;
  int minsize;
  double minlossimprove;

  NodeWorker(std::vector<Node>& in, std::vector<std::vector<Node>>& out,
             int minsize_, double minlossimprove_) :
    nodes_in(in), nodes_out(out), minsize(minsize_), minlossimprove(minlossimprove_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (size_t i = begin; i < end; i++) {
      if (nodes_in[i].leaf) nodes_out[i] = {nodes_in[i]};
      else nodes_out[i] = getchild(nodes_in[i], minsize, minlossimprove);
    }
  }
};

// ---------------- Top-level function ----------------
// [[Rcpp::export]]
List gettree_cpp(arma::vec y, arma::mat X, double minlossimprove=0.01, int minsize=1) {
  Node root = {y, X, false};
  std::vector<Node> current_gen = {root};
  std::vector<Node> next_gen;

  while (true) {
    // collect non-leaf nodes
    std::vector<Node> non_leaf_nodes;
    std::vector<size_t> indices;
    for (size_t i=0; i<current_gen.size(); i++) {
      if (!current_gen[i].leaf) {
        non_leaf_nodes.push_back(current_gen[i]);
        indices.push_back(i);
      }
    }
    if (non_leaf_nodes.empty()) break;

    std::vector<std::vector<Node>> optimized_children(non_leaf_nodes.size());
    NodeWorker worker(non_leaf_nodes, optimized_children, minsize, minlossimprove);
    RcppParallel::parallelFor(0, non_leaf_nodes.size(), worker);

    // flatten next generation
    next_gen.clear();
    size_t j = 0;
    for (size_t i=0; i<current_gen.size(); i++) {
      if (current_gen[i].leaf) next_gen.push_back(current_gen[i]);
      else {
        for (auto &child : optimized_children[j]) next_gen.push_back(child);
        j++;
      }
    }
    if (next_gen.size() == current_gen.size()) break;

    current_gen = next_gen;
  }

  List result(current_gen.size());
  for (size_t i=0; i<current_gen.size(); i++) {
    result[i] = List::create(
      _["intensity"] = current_gen[i].intensity,
      _["coord"] = current_gen[i].coord,
      _["leaf"] = current_gen[i].leaf
    );
  }

  return result;
}
