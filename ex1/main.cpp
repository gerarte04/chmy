#include "lib/argparse.hpp"
#include "decomposition.hpp"
#include "matrix.hpp"
#include "reader.hpp"
#include "solver.hpp"
#include <iostream>

int main(int argc, char **argv)
{
    argparse::ArgumentParser parser("QrDecompositionSolver");

    parser.add_argument("-i", "--input")
        .help("file with input matrix")
        .required();
    parser.add_argument("-o", "--output")
        .help("output file")
        .default_value("");
    parser.add_argument("--qr-method")
        .help("qr decomposition method (cholesky/householder)")
        .required()
        .choices("cholesky", "householder");
    parser.add_argument("-n", "--n-tests")
        .help("tests count")
        .scan<'u', unsigned int>()
        .default_value(5);

    try {
        parser.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string output = parser.get<std::string>("--output");

    if (output.length()) {
        assert(freopen(output.c_str(), "w", stdout) != nullptr);
    }

    vector<vector<double>> A;
    read_mtx(A, parser.get<std::string>("--input"));
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

    std::string method = parser.get<std::string>("--qr-method");

    if (method == "cholesky") {
        qr_decomposition_cholesky(A, Q, R);
    } else {
        qr_decomposition_householder(A, Q, R);
    }

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
    unsigned int n_tests = parser.get<unsigned int>("--n-tests");

    for (unsigned int i = 0; i < n_tests; ++i) {
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

    mean_disc /= n_tests;
    mean_err /= n_tests;
    std::cout << "Mean discrepancy: " << mean_disc << std::endl;
    std::cout << "Mean solution error: " << mean_err << std::endl;

    return 0;
}
