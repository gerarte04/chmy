#include "decomposition.hpp"
#include "matrix.hpp"
#include "reader.hpp"
#include "solver.hpp"
#include <iostream>

constexpr int N_TESTS = 5;

int main()
{
    assert(freopen("output.txt", "w", stdout) != nullptr);

    vector<vector<double>> A;
    read_mtx(A, "SLAU_var_2.csv");
    int n = A.size();

    vector<vector<double>> Q(n, vector<double>(n));
    vector<vector<double>> R(n, vector<double>(n));

    // for (auto &r : A) {
    //     for (auto v : r) {
    //         std::cout << v << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

#ifdef HOUSEHOLDER
    qr_decomposition_householder(A, Q, R);
#else
    qr_decomposition_holetskii(A, Q, R);
#endif
    for (auto &r : Q) {
        for (auto v : r) {
            std::cout << v << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (auto &r : R) {
        for (auto v : r) {
            std::cout << v << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Matrix difference norm: " << max_norm(A - Q * R) << std::endl << std::endl;

    double mean_disc{}, mean_err{};

    for (int i = 0; i < N_TESTS; ++i) {
        std::cout << "**Test " << i << "**\n";
        vector<double> x = generate_random(n);
        vector<double> f = A * x;
        vector<double> x_ = solve(R, solve(Q, f));

        double disc = max_norm(f - A * x_);
        double err = max_norm(x - x_);
        mean_disc += disc;
        mean_err += err;

        std::cout << "Max-norm of discrepancy: " << disc << std::endl;
        std::cout << "Max-norm of solution error: " << err << std::endl << std::endl;
    }

    mean_disc /= N_TESTS;
    mean_err /= N_TESTS;
    std::cout << "Mean discrepancy: " << mean_disc << std::endl;
    std::cout << "Mean solution error: " << mean_err << std::endl;

    return 0;
}
