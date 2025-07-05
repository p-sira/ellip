#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rg.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <functional>

// Helper to split a string by comma
std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream iss(line);
    while (std::getline(iss, token, ',')) {
        tokens.push_back(token);
    }
    return tokens;
}

// Generic process function
void process_carlson(
    const std::string& infile,
    const std::string& outfile,
    int nargs,
    std::function<double(const std::vector<double>&)> func
) {
    std::ifstream fin(infile);
    std::ofstream fout(outfile);
    std::string line;
    
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        auto tokens = split(line);
        if ((int)tokens.size() < nargs) continue;
        std::vector<double> args;
        for (int i = 0; i < nargs; ++i) {
            args.push_back(std::stod(tokens[i]));
        }
        double result = func(args);
        fout << std::setprecision(17);
        for (int i = 0; i < nargs; ++i) {
            fout << args[i];
            if (i < nargs - 1) fout << ",";
        }
        fout << "," << result << "\n";
    }
    std::cout << "Generated " << outfile << "." << std::endl;
}

int main() {
    using namespace boost::math;
    process_carlson(
        "../wolfram/elliprf_data.csv", "elliprf_data.csv", 3,
        [](const std::vector<double>& v) { return ellint_rf(v[0], v[1], v[2]); }
    );
    process_carlson(
        "../wolfram/elliprg_data.csv", "elliprg_data.csv", 3,
        [](const std::vector<double>& v) { return ellint_rg(v[0], v[1], v[2]); }
    );
    process_carlson(
        "../wolfram/elliprj_data.csv", "elliprj_data.csv", 4,
        [](const std::vector<double>& v) { return ellint_rj(v[0], v[1], v[2], v[3]); }
    );
    process_carlson(
        "../wolfram/elliprj_pv.csv", "elliprj_pv.csv", 4,
        [](const std::vector<double>& v) { return ellint_rj(v[0], v[1], v[2], v[3]); }
    );
    return 0;
}



