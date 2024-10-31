#pragma once

#include <assert.h>
#include <stdlib.h>
#include "matrix.hpp"

using std::vector;

vector<double> generate_random(size_t n)
{
    vector<double> x(n);

    for (size_t i = 0; i < n; ++i) {
        x[i] = (double) rand() / RAND_MAX * 2 - 1;
    }

    return x;
}

vector<double> solve(
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
