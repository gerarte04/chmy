#pragma once

#include <assert.h>
#include <random>
#include "matrix.hpp"

using std::vector;

std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());

vector<double> generate_random(size_t n)
{
    vector<double> x(n);

    for (size_t i = 0; i < n; ++i) {
        double x_norm = ((double) gen() - gen.min()) / (gen.max() - gen.min());
        x[i] = x_norm * 2 - 1;
    }

    return x;
}

vector<double> solve_lower_tr(
    const vector<vector<double>> &a,
    const vector<double> &f
)
{
    assert(
        a.size() && a[0].size() &&
        a.size() == a[0].size() && a[0].size() == f.size()
    );

    vector<double> x(f.size());

    for (size_t i = 0; i < f.size(); ++i) {
        double r = f[i];

        for (size_t j = 0; j < i; ++j) {
            r -= a[i][j] * x[j];
        }

        r /= a[i][i];
        x[i] = r;
    }

    return x;
}

vector<double> solve_upper_tr(
    const vector<vector<double>> &a,
    const vector<double> &f
)
{
    assert(
        a.size() && a[0].size() &&
        a.size() == a[0].size() && a[0].size() == f.size()
    );

    vector<double> x(f.size());

    for (int z = f.size() - 1; z >= 0; --z) {
        size_t i = z;
        double r = f[i];

        for (size_t j = i + 1; j < a[0].size(); ++j) {
            r -= a[i][j] * x[j];
        }

        r /= a[i][i];
        x[i] = r;
    }

    return x;
}

vector<double> solve_cholesky(
    const vector<vector<double>> &q,
    const vector<vector<double>> &r,
    const vector<double> &f
)
{
    return solve_upper_tr(r, solve_lower_tr(q, f));
}

vector<double> solve_householder(
    const vector<vector<double>> &q,
    const vector<vector<double>> &r,
    const vector<double> &f
)
{
    vector<vector<double>> qt(q.size(), vector<double>(q[0].size()));
    fill_transpose(q, qt);

    return solve_upper_tr(r, qt * f);
}
