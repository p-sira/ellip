/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#include <boost/math/special_functions/ellint_rj.hpp>
#include <math.h>
#include <iomanip>
#include <iostream>

using namespace boost::math;
int main()
{
    double ans = ellint_rj(2.22e-14, 2.22e-14, 1.6000000000000045, -2.22e-14);
    double rtol = std::abs((ans - -12.0312499999998) / ans / 2.22e-22);
    std::cout << std::setprecision(16) << ans << std::endl;
    std::cout << std::setprecision(16) << rtol << std::endl;

    return 0;
}