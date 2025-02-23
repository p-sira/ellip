/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

// Generate expected results for ellippi in double precision (f64),
// as the default Boost test data uses double promotion to long double (f128).

#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/ellip3.hpp>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "util.cpp"

using namespace std;
using namespace boost::math::policies;

int main()
{
    ifstream f_in("../tests/data/boost/ellippi2_data.txt");
    ofstream f_out("../tests/data/boost/ellippi2_data_f64.txt");
    if (!f_in.is_open() || !f_out.is_open())
    {
        cerr << "Cannot open file" << endl;
        return 1;
    }

    typedef policy<promote_double<false>> f64_policy;

    string line;
    while (getline(f_in, line))
    {
        vector inp = split_words(line);
        double k = stod(inp[1]);
        double v = stod(inp[0]);
        double ans = boost::math::ellip3(k, v, f64_policy());
        f_out << line << "    " << setprecision(17) << ans << "\n";
    }
    return 0;
}
