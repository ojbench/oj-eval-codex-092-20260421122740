// Problem 092 - Resistive Network implementation
#ifndef SRC_HPP
#define SRC_HPP

#include <vector>
#include <algorithm>
#include <iostream>
#include "fraction.hpp"

class resistive_network {
private:
    int n; // nodes
    int m; // edges
    struct Edge { int u, v; fraction r; };
    std::vector<Edge> edges;

    // Build Laplacian matrix L with conductances g = 1/r
    void build_laplacian(std::vector<std::vector<fraction>> &L) const {
        L.assign(n, std::vector<fraction>(n, fraction(0)));
        for (const auto &e : edges) {
            int a = e.u - 1;
            int b = e.v - 1;
            fraction g = fraction(1) / e.r;
            L[a][a] = L[a][a] + g;
            L[b][b] = L[b][b] + g;
            L[a][b] = L[a][b] - g;
            L[b][a] = L[b][a] - g;
        }
    }

    // Solve linear system A x = b using fraction Gaussian elimination
    static std::vector<fraction> solve_linear(std::vector<std::vector<fraction>> A,
                                              std::vector<fraction> b) {
        int n = (int)A.size();
        if (n == 0) throw matrix_error();
        for (auto &row : A) if ((int)row.size() != n) throw matrix_error();
        if ((int)b.size() != n) throw matrix_error();

        // Augment matrix [A|b]
        for (int i = 0; i < n; ++i) A[i].push_back(b[i]);

        for (int col = 0, row = 0; col < n && row < n; ++col) {
            int sel = -1;
            for (int i = row; i < n; ++i) {
                if (!(A[i][col] == fraction(0))) { sel = i; break; }
            }
            if (sel == -1) throw matrix_error(); // singular
            if (sel != row) std::swap(A[sel], A[row]);

            fraction piv = A[row][col];
            for (int j = col; j <= n; ++j) A[row][j] = A[row][j] / piv;

            for (int i = 0; i < n; ++i) {
                if (i == row) continue;
                fraction factor = A[i][col];
                if (factor == fraction(0)) continue;
                for (int j = col; j <= n; ++j) {
                    A[i][j] = A[i][j] - factor * A[row][j];
                }
            }
            ++row;
        }

        std::vector<fraction> x(n, fraction(0));
        for (int i = 0; i < n; ++i) x[i] = A[i][n];
        return x;
    }

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        n = interface_size_;
        m = connection_size_;
        edges.reserve(m);
        for (int i = 0; i < m; ++i) {
            edges.push_back(Edge{from[i], to[i], resistance[i]});
        }
    }

    // Equivalent resistance between two nodes with 1A injection and ground at node n
    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int a = interface_id1;
        int b = interface_id2;
        if (a == b) return fraction(0);
        std::vector<std::vector<fraction>> L;
        build_laplacian(L);
        int dim = n - 1;
        std::vector<std::vector<fraction>> A(dim, std::vector<fraction>(dim, fraction(0)));
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                A[i][j] = L[i][j];
        std::vector<fraction> I(dim, fraction(0));
        if (a <= dim) I[a - 1] = I[a - 1] + fraction(1);
        if (b <= dim) I[b - 1] = I[b - 1] - fraction(1);
        std::vector<fraction> U = solve_linear(A, I);
        fraction ua = (a == n) ? fraction(0) : U[a - 1];
        fraction ub = (b == n) ? fraction(0) : U[b - 1];
        return ua - ub;
    }

    // Voltage at node id given currents (u_n = 0)
    fraction get_voltage(int id, fraction current[]) {
        std::vector<std::vector<fraction>> L;
        build_laplacian(L);
        int dim = n - 1;
        std::vector<std::vector<fraction>> A(dim, std::vector<fraction>(dim, fraction(0)));
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                A[i][j] = L[i][j];
        std::vector<fraction> I(dim, fraction(0));
        for (int i = 0; i < dim; ++i) I[i] = current[i];
        std::vector<fraction> U = solve_linear(A, I);
        if (id == n) return fraction(0);
        return U[id - 1];
    }

    // Network power given node voltages: u^T L u
    fraction get_power(fraction voltage[]) {
        std::vector<std::vector<fraction>> L;
        build_laplacian(L);
        std::vector<fraction> u(n);
        for (int i = 0; i < n; ++i) u[i] = voltage[i];
        std::vector<fraction> Lu(n, fraction(0));
        for (int i = 0; i < n; ++i) {
            fraction s(0);
            for (int j = 0; j < n; ++j) s = s + L[i][j] * u[j];
            Lu[i] = s;
        }
        fraction utLu(0);
        for (int i = 0; i < n; ++i) utLu = utLu + u[i] * Lu[i];
        return utLu;
    }
};

#endif // SRC_HPP

