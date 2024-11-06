#pragma once

#include <vector>
#include "lib/csv.hpp"

using namespace csv;
using std::vector;

void read_mtx(vector<vector<double>> &mtx, std::string &&input_file)
{
    CSVReader reader(input_file, CSVFormat().no_header());

    for (CSVRow &row : reader) {
        vector<double> r;

        for (CSVField &field : row) {
            r.push_back(field.get<double>());
        }

        mtx.push_back(r);
    }
}
