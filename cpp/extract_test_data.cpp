/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#include <iostream>
#include <fstream>
#include <string>
#include <regex>

// Extract txt file from Boost Math test data ipp
void extract_ipp_data(const std::string &input_file, const std::string &output_file)
{
    std::ifstream in(input_file);
    std::ofstream out(output_file);
    std::string line;

    while (std::getline(in, line))
    {
        if (line.find("static const std::array") != std::string::npos)
        {
            break;
        }
    }

    std::regex number_pattern(R"(SC_\(([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\))");

    while (std::getline(in, line))
    {
        if (line.find("}}") != std::string::npos && line.length() < 5)
        {
            break;
        }

        if (line.find("SC_") != std::string::npos)
        {
            std::string cleaned_line;
            auto numbers_begin = std::sregex_iterator(line.begin(), line.end(), number_pattern);
            auto numbers_end = std::sregex_iterator();

            for (auto i = numbers_begin; i != numbers_end; ++i)
            {
                std::smatch match = *i;
                cleaned_line += match[1].str() + " ";
            }

            if (!cleaned_line.empty())
            {
                out << cleaned_line.substr(0, cleaned_line.length() - 1) << "\n";
            }
        }
    }
}

int main()
{
    extract_ipp_data("ellint_d2_data.ipp", "../tests/data/boost/ellipdinc_data.txt");
    return 0;
}
