#pragma once

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <vector>

using std::vector;

double l2_norm(const vector<double> &v)
{
    assert(v.size());
    double s{};

    for (size_t i = 0; i < v.size(); ++i) {
        s += pow(v[i], 2);
    }

    return sqrt(s);
}

double max_norm(const vector<double> &v)
{
    assert(v.size());
    return *std::max_element(v.begin(), v.end());
}

double max_norm(const vector<vector<double>> &mtx)
{
    assert(mtx.size() && mtx[0].size());

    auto sum_row = [](const vector<double> &row) {
        double sum{};

        for (const double &a : row) {
            sum += a;
        }

        return sum;
    };

    double max_row = sum_row(mtx.front());

    for (const auto &row : mtx) {
        max_row = std::max(max_row, sum_row(row));
    }

    return max_row;
}

vector<double> operator-(const vector<double> &a, const vector<double> &b) 
{
    assert(a.size() && a.size() == b.size());

    vector<double> new_a(a.begin(), a.end());
    
    for (size_t i = 0; i < a.size(); ++i) {
        new_a[i] -= b[i];
    }

    return new_a;
}

vector<vector<double>> operator-(
    const vector<vector<double>> &a,
    const vector<vector<double>> &b
)
{
    assert(
        a.size() == b.size() && a.size() && b.size() &&
        a[0].size() && a[0].size() == b[0].size()
    );

    vector<vector<double>> new_a(a.begin(), a.end());

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a[0].size(); ++j) {
            new_a[i][j] -= b[i][j];
        }
    }

    return new_a;
}

vector<vector<double>> operator*(
    const vector<vector<double>> &a,
    const vector<vector<double>> &b
)
{
    assert(
        a.size() && b.size() && a[0].size() &&
        a[0].size() == b.size()
    );

    vector<vector<double>> c(a.size(), vector<double>(b[0].size()));

    for (size_t i = 0; i < c.size(); ++i) {
        for (size_t j = 0; j < c[0].size(); ++j) {
            double sum{};

            for (size_t k = 0; k < a[0].size(); ++k) {
                sum += a[i][k] * b[k][j];
            }

            c[i][j] = sum;
        }
    }

    return c;
}

void fill_transpose(
    const vector<vector<double>> &a,
    vector<vector<double>> &at
)
{
    assert(
        a.size() && a[0].size() && at.size() &&
        a.size() == at[0].size() && a[0].size() == at.size()
    );

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a[0].size(); ++j) {
            at[j][i] = a[i][j];
        }
    }
}

vector<vector<double>>& operator*(
    vector<vector<double>> &&a,
    double &&alpha    
)
{
    assert(a.size() && a[0].size());

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a[0].size(); ++j) {
            a[i][j] *= alpha;
        }
    }

    return a;
}

vector<vector<double>> uut_mul(vector<double> &u)
{
    assert(u.size());
    vector<vector<double>> a(u.size(), vector<double>(u.size()));

    for (size_t i = 0; i < u.size(); ++i) {
        for (size_t j = 0; j < u.size(); ++j) {
            a[i][j] = u[i] * u[j];
        }
    }

    return a;
}

vector<vector<double>> one_mtx(size_t n)
{
    vector<vector<double>> a(n, vector<double>(n));

    for (size_t i = 0; i < n; ++i) {
        a[i][i] = 1;
    }

    return a;
}

vector<double> operator*(
    const vector<vector<double>> &a,
    const vector<double> &x
)
{
    assert(a.size() && a[0].size() == x.size());

    vector<double> y(x.size());

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a[0].size(); ++j) {
            y[i] += a[i][j] * x[j];
        }
    }

    return y;
}
