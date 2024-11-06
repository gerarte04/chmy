#pragma once

#include <assert.h>
#include <math.h>
#include <vector>
#include "matrix.hpp"

using std::vector;

void qr_decomposition_cholesky(
    vector<vector<double>> &a,
    vector<vector<double>> &q,
    vector<vector<double>> &r
)
{
    assert(
        a.size() && a[0].size() && q.size() && r.size() &&
        a.size() == q.size() && a[0].size() == q[0].size() &&
        a.size() == r.size() && a[0].size() == r[0].size()
    );

    q[0][0] = sqrt(a[0][0]);

    for (size_t j = 1; j < a.size(); ++j) {
        q[j][0] = a[j][0] / q[0][0];
    }

    for (size_t i = 1; i < a[0].size(); ++i) {
        double s1{};

        for (size_t j = 0; j < i; ++j) {
            s1 += pow(q[i][j], 2);
        }

        q[i][i] = sqrt(a[i][i] - s1);

        for (size_t j = i + 1; j < a.size(); ++j) {
            double s2{};

            for (size_t p = 0; p < i; ++p) {
                s2 += q[i][p] * q[j][p];
            }

            q[j][i] = (a[j][i] - s2) / q[i][i];
        }
    }

    fill_transpose(q, r);
}

void qr_decomposition_householder(
    vector<vector<double>> &ca,
    vector<vector<double>> &q,
    vector<vector<double>> &r
)
{
    assert(
        ca.size() && ca[0].size() && q.size() && r.size() &&
        ca.size() == q.size() && ca[0].size() == q[0].size() &&
        ca.size() == r.size() && ca[0].size() == r[0].size()
    );

    int n = ca.size();
    vector<vector<double>> a(ca.begin(), ca.end());

    for (int k = 0; k < n - 1; ++k) {
        vector<double> x(n);

        for (int i = 0; i < n; ++i) {
            x[i] = (i < k) ? 0 : a[i][k];
        }

        vector<double> e(n);
        e[k] = l2_norm(x);
        vector<double> u = x - e;

        vector<vector<double>> p = one_mtx(n) - uut_mul(u) * (2 / pow(l2_norm(u), 2));
        q = q * p;
        a = p * a;
        r = a;
    }
}
