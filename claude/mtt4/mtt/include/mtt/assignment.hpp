#pragma once
#include "types.hpp"
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cassert>

namespace mtt {

// ─────────────────────────────────────────────────────────────────────────────
// 2-D Hungarian Algorithm (Kuhn-Munkres)
// Solves the minimum-cost assignment of M measurements to N tracks.
//
// cost[i][j] = assignment cost for measurement i → track j.
// Returns: assignment[i] = j  (or -1 if measurement i is unassigned).
// Complexity: O(min(M,N)²·max(M,N))
// ─────────────────────────────────────────────────────────────────────────────
class HungarianSolver {
public:
    using CostMatrix = std::vector<std::vector<double>>;

    // Returns assignments vector of length n_meas.
    // assignment[i] = track_index, or -1 if unassigned.
    static std::vector<int> solve(const CostMatrix& cost,
                                  double unassigned_cost = 1e9) {
        const int n_meas  = static_cast<int>(cost.size());
        if (n_meas == 0) return {};
        const int n_track = static_cast<int>(cost[0].size());

        // Square up the cost matrix with dummy entries
        const int N = std::max(n_meas, n_track);
        std::vector<std::vector<double>> C(N, std::vector<double>(N, unassigned_cost));
        for (int i = 0; i < n_meas;  ++i)
        for (int j = 0; j < n_track; ++j)
            C[i][j] = cost[i][j];

        // Munkres
        std::vector<int> u(N+1), v(N+1), p(N+1), way(N+1);
        std::fill(u.begin(), u.end(), 0);
        std::fill(v.begin(), v.end(), 0);
        std::fill(p.begin(), p.end(), 0);

        for (int i = 1; i <= N; ++i) {
            p[0] = i;
            int j0 = 0;
            std::vector<double> minv(N+1, std::numeric_limits<double>::max());
            std::vector<bool> used(N+1, false);
            do {
                used[j0] = true;
                int i0 = p[j0], j1 = -1;
                double delta = std::numeric_limits<double>::max();
                for (int j = 1; j <= N; ++j) {
                    if (!used[j]) {
                        double cur = C[i0-1][j-1] - u[i0] - v[j];
                        if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                        if (minv[j] < delta) { delta = minv[j]; j1 = j; }
                    }
                }
                for (int j = 0; j <= N; ++j) {
                    if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                    else minv[j] -= delta;
                }
                j0 = j1;
            } while (p[j0] != 0);
            do {
                int j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
            } while (j0);
        }

        std::vector<int> assign(n_meas, -1);
        for (int j = 1; j <= N; ++j) {
            if (p[j] >= 1 && p[j] <= n_meas && j-1 < n_track)
                assign[p[j]-1] = j-1;
        }
        return assign;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Multi-dimensional (S-D) Assignment via Lagrangian Relaxation.
//
// Fuses measurement-to-track associations from S ≤ 4 sensors simultaneously.
// The NP-hard S-D assignment is relaxed into a sequence of 2-D Hungarian
// problems. Lagrange multipliers penalise constraint violations and are
// updated iteratively until convergence.
//
// Reference: Poore & Rijavec, IEEE TAC 1993; Deb et al., PAMI 1997.
// ─────────────────────────────────────────────────────────────────────────────
struct AssignmentProblem {
    // cost_matrix[sensor_idx][meas_idx][track_idx] = Mahalanobis² cost
    std::vector<std::vector<std::vector<double>>> cost_matrices;
    int n_tracks;
    int n_sensors;

    // unassigned penalty (> gate threshold to disfavour clutter associations)
    double unassigned_penalty = 20.0;
};

struct AssignmentResult {
    // For each sensor: assignment[meas_idx] = track_idx (or -1)
    std::vector<std::vector<int>> sensor_assignments;
    double total_cost;
};

class LagrangianMDA {
public:
    static constexpr int MAX_ITER       = 50;
    static constexpr double STEP_INIT   = 1.0;
    static constexpr double STEP_DECAY  = 0.95;
    static constexpr double TOL         = 1e-4;

    static AssignmentResult solve(const AssignmentProblem& prob) {
        const int S = prob.n_sensors;
        const int T = prob.n_tracks;

        AssignmentResult result;
        result.sensor_assignments.resize(S);

        if (S == 0 || T == 0) { result.total_cost = 0.0; return result; }

        // ── Single-sensor: direct Hungarian ───────────────────────────────────
        if (S == 1) {
            result.sensor_assignments[0] =
                HungarianSolver::solve(prob.cost_matrices[0],
                                       prob.unassigned_penalty);
            result.total_cost = 0.0;
            return result;
        }

        // ── Multi-sensor: Lagrangian relaxation ───────────────────────────────
        // Lagrange multipliers λ[t] for each track (one per track constraint).
        std::vector<double> lambda(T, 0.0);

        // Best feasible solution bookkeeping
        AssignmentResult best_result;
        best_result.sensor_assignments.resize(S);
        double best_cost = std::numeric_limits<double>::max();

        double step = STEP_INIT;

        for (int iter = 0; iter < MAX_ITER; ++iter) {
            // ── Solve S independent 2-D sub-problems with perturbed costs ─────
            std::vector<std::vector<int>> sub_assign(S);
            double lower_bound = 0.0;

            for (int s = 0; s < S; ++s) {
                // Perturb cost by Lagrange multipliers
                auto C = prob.cost_matrices[s];   // copy
                for (auto& row : C) {
                    for (int t = 0; t < T && t < static_cast<int>(row.size()); ++t)
                        row[t] += lambda[t];
                }
                sub_assign[s] = HungarianSolver::solve(C, prob.unassigned_penalty);
                // Accumulate lower bound
                for (std::size_t m = 0; m < sub_assign[s].size(); ++m) {
                    int t = sub_assign[s][m];
                    if (t >= 0 && t < T)
                        lower_bound += prob.cost_matrices[s][m][t];
                    else
                        lower_bound += prob.unassigned_penalty;
                }
            }
            lower_bound -= std::accumulate(lambda.begin(), lambda.end(), 0.0)
                           * static_cast<double>(S - 1);

            // ── Compute sub-gradient ──────────────────────────────────────────
            std::vector<double> grad(T, -static_cast<double>(S));
            for (int s = 0; s < S; ++s) {
                for (int t : sub_assign[s]) {
                    if (t >= 0 && t < T) grad[t] += 1.0;
                }
            }

            // ── Update λ (sub-gradient ascent to maximise lower bound) ─────────
            double norm2 = 0.0;
            for (double g : grad) norm2 += g*g;
            if (norm2 < TOL) break;

            for (int t = 0; t < T; ++t)
                lambda[t] = std::max(0.0, lambda[t] + step * grad[t]);

            step *= STEP_DECAY;

            // ── Feasibility check (each track assigned at most once) ───────────
            bool feasible = true;
            std::vector<int> track_count(T, 0);
            for (int s = 0; s < S; ++s) {
                for (int t : sub_assign[s]) {
                    if (t >= 0 && t < T) {
                        ++track_count[t];
                        if (track_count[t] > S) { feasible = false; }
                    }
                }
            }

            if (feasible && lower_bound < best_cost) {
                best_cost = lower_bound;
                best_result.sensor_assignments = sub_assign;
                best_result.total_cost = lower_bound;
            }
        }

        return best_result;
    }
};

}  // namespace mtt
