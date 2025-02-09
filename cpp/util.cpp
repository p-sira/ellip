/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#include <vector>
#include <string>
#include <iterator>
#include <sstream>

std::vector<std::string> split_words(const std::string &input)
{
    std::istringstream buffer(input);
    return {std::istream_iterator<std::string>(buffer), std::istream_iterator<std::string>()};
}