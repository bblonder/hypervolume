// SUGGESTED TAB SIZE: 4
#ifndef _KDTREE_H_
#define _KDTREE_H_

#include <vector>
#include <cstddef>
#include <math.h>
#include <utility>
#include "MyHeaps.h"
#include "float.h"

using namespace std;

/// The root node is stored in position 0 of nodes
#define ROOT 0

class Node {
public:
  double  key;   // the key (value along k-th dimension) of the split
  int     LIdx; // the index to the left cell  (-1 if none)
  int     RIdx; // the index to the right cell (-1 if none)

  // A pointer back to the structure data of input points,
  // but ONLY if the node is a LEAF, otherwise value is (-1).
  int     pIdx;

  inline bool isLeaf() const { return pIdx >= 0; }

  Node() : key(-1), LIdx(-1), RIdx(-1), pIdx(-1) {}
};

class KDTree {
  // Point coordinates stored ROW-MAJOR in a single contiguous buffer:
  //   coord(i, j) == points[(size_t) i * ndim + j]
  // Previously vector<vector<double>>, one heap allocation per point,
  // scattering the tree's working set across the heap.
private: vector<double> points;

  // Tree nodes stored BY VALUE in one contiguous buffer. Previously
  // vector<Node*> with `new Node()` per split, ~2n heap allocations for
  // an n-point tree and a cache miss per traversal step.
private: vector<Node> nodes;

private: int ndim;                // Data dimensionality
private: int npoints;             // Number of points
private: vector<int> workarray;   // Used in tree construction

  // Data used ONLY for "kNN" search (k-nearest neighbors)
private: vector<double> Bmin;              // bounding box lower bound
private: vector<double> Bmax;              // bounding box upper bound
private: MaxHeap<double> pq;               // <key,idx> = <distance, node idx>
private: int k;                            // number of records to search for
private: bool terminate_search;            // true if k points have been found

  // ---- point access helpers -------------------------------------------------

private:
  inline const double* point(int i) const {
    return &points[(size_t) i * (size_t) ndim];
  }
  inline double coord(int i, int j) const {
    return points[(size_t) i * (size_t) ndim + (size_t) j];
  }
  // Squared L2 distance between query row q (length ndim) and stored point i.
  inline double distance_squared_to_point(const double* q, int i) const {
    const double* pi = point(i);
    double d = 0;
    for (int j = 0; j < ndim; j++) {
      double diff = q[j] - pi[j];
      d += diff * diff;
    }
    return d;
  }
  // Squared L2 distance between two arbitrary rows (both length ndim).
  inline double distance_squared(const double* a, const double* b) const {
    double d = 0;
    for (int j = 0; j < ndim; j++) {
      double diff = a[j] - b[j];
      d += diff * diff;
    }
    return d;
  }
  inline bool point_lies_in_range(int i, const vector<double>& pMin,
                                  const vector<double>& pMax) const {
    const double* p = point(i);
    for (int j = 0; j < ndim; j++) {
      if (p[j] < pMin[j] || p[j] > pMax[j]) return false;
    }
    return true;
  }

public:
  /// Nothing to release explicitly; vector<Node> owns the storage.
  ~KDTree() {}

  /// Number of points in the kd-tree.
  inline int size() const { return npoints; }
  /// Dimensionality of the points in the kd-tree.
  inline int ndims() const { return ndim; }

  /**
   * Creates a KDtree filled with the provided data. The buffer is
   * ROW-MAJOR: element (row i, col j) at flat_row_major[i*ncol + j].
   * Ownership of the buffer is taken (moved in) so the caller need not
   * keep it alive.
   */
  KDTree(vector<double>&& flat_row_major, int npoints_, int ndim_)
    : points(std::move(flat_row_major)),
      ndim(ndim_),
      npoints(npoints_)
  {
    nodes.reserve((size_t) 2 * npoints);  // 2n-1 nodes upper bound
    workarray.assign(npoints, -1);        // -1 invariant per build_recursively

    // Build sorted-index tables in each dimension via MinHeap.
    vector< MinHeap<double> > heaps(ndim, npoints);
    for (int d = 0; d < ndim; d++) {
      for (int p = 0; p < npoints; p++) {
        heaps[d].push(coord(p, d), p);
      }
    }
    // indexes[dim][i]: in dimension dim, index of the i-th smallest element.
    vector< vector<int> > indexes(ndim, vector<int>(npoints, 0));
    for (int d = 0; d < ndim; d++) {
      heaps[d].heapsort(indexes[d]);
    }
    // srtidx[dim][j]: in dimension dim, ordering number of point j.
    vector< vector<int> > srtidx(ndim, vector<int>(npoints, 0));
    for (int d = 0; d < ndim; d++) {
      for (int i = 0; i < npoints; i++) {
        srtidx[d][indexes[d][i]] = i;
      }
    }

    // First partition is on every single point.
    vector<int> pidx(npoints, 0);
    for (int i = 0; i < npoints; i++) pidx[i] = i;

    build_recursively(srtidx, pidx, 0);
  }

  /**
   * Heapsort algorithm used by fast KDtree construction.
   * Copied almost verbatim from Wikipedia's Heapsort page.
   */
private:
  int heapsort(int dim, vector<int>& idx, int len) {
    unsigned int n = len;
    unsigned int i = len / 2;
    unsigned int parent, child;
    int t;

    for (;;) {
      if (i > 0) {
        i--;
        t = idx[i];
      } else {
        n--;
        if (n == 0) return 0;
        t = idx.at(n);
        idx[n] = idx[0];
      }

      parent = i;
      child = i * 2 + 1;

      while (child < n) {
        if ((child + 1 < n) && (coord(idx[child + 1], dim) > coord(idx[child], dim))) {
          child++;
        }
        if (coord(idx[child], dim) > coord(t, dim)) {
          idx[parent] = idx[child];
          parent = child;
          child = parent * 2 + 1;
        } else {
          break;
        }
      }
      idx[parent] = t;
    }
    return 0;
  }

  /**
   * Recursively performs median splits along dimension "dim" using the
   * pre-prepared per-dimension sort orderings. Memory-friendly version.
   */
private:
  int build_recursively(vector< vector<int> >& sortidx, vector<int>& pidx, int dim) {
    // Stop condition: one point left => this is a leaf.
    if (pidx.size() == 1) {
      int nodeIdx = (int) nodes.size();
      nodes.emplace_back();
      Node& node = nodes.back();
      node.LIdx = -1;
      node.RIdx = -1;
      node.pIdx = pidx[0];
      node.key  = 0;
      return nodeIdx;
    }

    // allocate the child index vectors
    vector<int> Larray;
    vector<int> Rarray;
    Larray.reserve(pidx.size() / 2 + (pidx.size() % 2 == 0 ? 0 : 1));
    Rarray.reserve(pidx.size() / 2);

    // Invariant: workarray is all -1 on entry. The heapsort branch does
    // not read workarray at all, so we skip the fill entirely there.
    // The scan branch fills only the pidx.size() positions it needs and
    // restores them at the end, so the invariant holds across sibling
    // recursive calls. (Proposal #1 in the perf review notes.)
    int pivot = -1;
    if ((double) pidx.size() * log((double) pidx.size()) < npoints) {
      Larray.resize(pidx.size() / 2 + (pidx.size() % 2 == 0 ? 0 : 1), -1);
      Rarray.resize(pidx.size() / 2, -1);
      heapsort(dim, pidx, (int) pidx.size());
      std::copy(pidx.begin(), pidx.begin() + Larray.size(), Larray.begin());
      std::copy(pidx.begin() + Larray.size(), pidx.end(),   Rarray.begin());
      pivot = pidx[(pidx.size() - 1) / 2];
    } else {
      for (size_t i = 0; i < pidx.size(); i++) {
        workarray[sortidx[dim][pidx[i]]] = pidx[i];
      }

      unsigned int TH = (unsigned int) (pidx.size() - 1) / 2;
      unsigned int cnt = 0;
      for (int i = 0; i < npoints; i++) {
        if (workarray[i] == -1) continue;
        if (cnt == TH) {
          pivot = workarray[i];
          Larray.push_back(workarray[i]);
        } else if (cnt > TH) {
          Rarray.push_back(workarray[i]);
        } else {
          Larray.push_back(workarray[i]);
        }
        cnt++;
        if (cnt > pidx.size()) break;
      }

      // Restore the -1 invariant on the positions we wrote.
      for (size_t i = 0; i < pidx.size(); i++) {
        workarray[sortidx[dim][pidx[i]]] = -1;
      }
    }

    // CREATE THE INTERNAL NODE. Reserve its slot BEFORE recursing so its
    // index is stable while children are added; then set fields once we
    // know LIdx and RIdx. (nodes may reallocate as children are added,
    // so we cannot hold a Node& across recursive calls.)
    int nodeIdx = (int) nodes.size();
    nodes.emplace_back();
    // Note: we read points before the recursive calls because we know
    // pivot is a valid index into points that won't be invalidated.
    double node_key = coord(pivot, dim);

    int lchild = build_recursively(sortidx, Larray, (dim + 1) % ndim);
    int rchild = build_recursively(sortidx, Rarray, (dim + 1) % ndim);

    // Re-fetch node reference: nodes may have reallocated during recursion.
    Node& node = nodes[nodeIdx];
    node.pIdx = -1;      // not a leaf
    node.key  = node_key;
    node.LIdx = lchild;
    node.RIdx = rchild;

    return nodeIdx;
  }

public:
  /**
   * k-NN query: computes the k closest points in the database to a given
   * point and returns their indexes and distances.
   */
  void k_closest_points(const vector<double>& Xq, int k_,
                        vector<int>& idxs, vector<double>& distances) {
    Bmin = vector<double>(ndim, -DBL_MAX);
    Bmax = vector<double>(ndim, +DBL_MAX);
    this->k = k_;
    this->terminate_search = false;

    knn_search(Xq);

    int N = (int) pq.size();
    for (int i = 0; i < N; i++) {
      pair<double, int> topel = pq.top();
      pq.pop();
      if (i >= N - k_) {
        idxs.push_back(topel.second);
        distances.push_back(sqrt(topel.first));
      }
    }
    // NOTE: the original had two consecutive std::reverse() calls here that
    // cancel out. Preserving the behaviour (a no-op) as an explicit no-op.
  }

private:
  void knn_search(const vector<double>& Xq, int nodeIdx = 0, int dim = 0) {
    Node& node = nodes[nodeIdx];
    double temp;

    if (node.isLeaf()) {
      double distance = distance_squared_to_point(Xq.data(), node.pIdx);

      if ((int) pq.size() == k && pq.top().first > distance) {
        pq.pop();
        pq.push(distance, node.pIdx);
      } else if ((int) pq.size() < k) {
        pq.push(distance, node.pIdx);
      }
      return;
    }

    if (Xq[dim] <= node.key) {
      temp = Bmax[dim]; Bmax[dim] = node.key;
      knn_search(Xq, node.LIdx, (dim + 1) % ndim);
      Bmax[dim] = temp;
    } else {
      temp = Bmin[dim]; Bmin[dim] = node.key;
      knn_search(Xq, node.RIdx, (dim + 1) % ndim);
      Bmin[dim] = temp;
    }
    if (Xq[dim] <= node.key) {
      temp = Bmin[dim]; Bmin[dim] = node.key;
      if (bounds_overlap_ball(Xq))
        knn_search(Xq, node.RIdx, (dim + 1) % ndim);
      Bmin[dim] = temp;
    } else {
      temp = Bmax[dim]; Bmax[dim] = node.key;
      if (bounds_overlap_ball(Xq))
        knn_search(Xq, node.LIdx, (dim + 1) % ndim);
      Bmax[dim] = temp;
    }
  }

  bool ball_within_bounds(const vector<double>& Xq) {
    double best_dist = sqrt(pq.top().first);
    for (int d = 0; d < ndim; d++) {
      if (fabs(Xq[d] - Bmin[d]) < best_dist || fabs(Xq[d] - Bmax[d]) < best_dist) {
        return false;
      }
    }
    return true;
  }

  double bounds_overlap_ball(const vector<double>& Xq) {
    if ((int) pq.size() < k) return true;

    double sum = 0;
    double best_dist_sq = pq.top().first;
    for (int d = 0; d < ndim; d++) {
      if (Xq[d] < Bmin[d]) {
        sum += (Xq[d] - Bmin[d]) * (Xq[d] - Bmin[d]);
        if (sum > best_dist_sq) return false;
      } else if (Xq[d] > Bmax[d]) {
        sum += (Xq[d] - Bmax[d]) * (Xq[d] - Bmax[d]);
        if (sum > best_dist_sq) return false;
      }
    }
    return true;
  }

public:
  /// Computes the closest point in the set to the query point.
  int closest_point(const vector<double>& p) {
    int dim = 0;
    int leafIdx = 0;
    for (;;) {
      const Node& leaf = nodes[leafIdx];
      if (leaf.pIdx >= 0) break;
      if (p[dim] <= leaf.key) leafIdx = leaf.LIdx;
      else                    leafIdx = leaf.RIdx;
      dim = (dim + 1) % ndim;
    }
    int leaf_pIdx = nodes[leafIdx].pIdx;
    double cdistsq = distance_squared_to_point(p.data(), leaf_pIdx);
    int closest_neighbor = leaf_pIdx;
    check_border_distance(ROOT, 0, p, cdistsq, closest_neighbor);
    return closest_neighbor;
  }

private:
  void check_border_distance(int nodeIdx, int dim, const vector<double>& pnt,
                             double& cdistsq, int& idx) {
    Node& node = nodes[nodeIdx];
    if (node.pIdx >= 0) {
      double dsq = distance_squared_to_point(pnt.data(), node.pIdx);
      if (dsq < cdistsq) {
        cdistsq = dsq;
        idx = node.pIdx;
      }
      return;
    }

    double ndistsq = (node.key - pnt[dim]) * (node.key - pnt[dim]);

    if (ndistsq > cdistsq) {
      if (node.key > pnt[dim])
        check_border_distance(node.LIdx, (dim + 1) % ndim, pnt, cdistsq, idx);
      else
        check_border_distance(node.RIdx, (dim + 1) % ndim, pnt, cdistsq, idx);
    } else {
      check_border_distance(node.LIdx, (dim + 1) % ndim, pnt, cdistsq, idx);
      check_border_distance(node.RIdx, (dim + 1) % ndim, pnt, cdistsq, idx);
    }
  }

public:
  /**
   * Query all points at distance less than or equal to radius from point.
   */
  void ball_query(const vector<double>& point_, const double radius,
                  vector<int>& idxsInRange, vector<double>& distances) {
    vector<double> pmin(ndim, 0);
    vector<double> pmax(ndim, 0);
    for (int dim = 0; dim < ndim; dim++) {
      pmin[dim] = point_[dim] - radius;
      pmax[dim] = point_[dim] + radius;
    }
    ball_bbox_query(ROOT, pmin, pmax, idxsInRange, distances, point_, radius * radius, 0);
  }

private:
  void ball_bbox_query(int nodeIdx, vector<double>& pmin, vector<double>& pmax,
                       vector<int>& inrange_idxs, vector<double>& distances,
                       const vector<double>& query, const double& radiusSquared,
                       int dim = 0) {
    Node& node = nodes[nodeIdx];

    if (node.isLeaf()) {
      double distance = distance_squared_to_point(query.data(), node.pIdx);
      if (distance <= radiusSquared) {
        inrange_idxs.push_back(node.pIdx);
        distances.push_back(sqrt(distance));
        return;
      }
    } else {
      if (node.key >= pmin[dim] && node.LIdx != -1)
        ball_bbox_query(node.LIdx, pmin, pmax, inrange_idxs, distances, query, radiusSquared, (dim + 1) % ndim);
      if (node.key <= pmax[dim] && node.RIdx != -1)
        ball_bbox_query(node.RIdx, pmin, pmax, inrange_idxs, distances, query, radiusSquared, (dim + 1) % ndim);
    }
  }

public:
  /**
   * k-dimensional Range query: returns indexes of points within the
   * bounding box [pmin, pmax].
   */
  void range_query(const vector<double>& pmin, const vector<double>& pmax,
                   vector<int>& inrange_idxs, int nodeIdx = 0, int dim = 0) {
    Node& node = nodes[nodeIdx];

    if (node.isLeaf()) {
      if (point_lies_in_range(node.pIdx, pmin, pmax)) {
        inrange_idxs.push_back(node.pIdx);
        return;
      }
    } else {
      if (node.key >= pmin[dim] && node.LIdx != -1)
        range_query(pmin, pmax, inrange_idxs, node.LIdx, (dim + 1) % ndim);
      if (node.key <= pmax[dim] && node.RIdx != -1)
        range_query(pmin, pmax, inrange_idxs, node.RIdx, (dim + 1) % ndim);
    }
  }
};

#endif
