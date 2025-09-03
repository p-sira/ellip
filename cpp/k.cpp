/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#include <boost/math/ccmath/sqrt.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <iostream>

using namespace boost::math;
using namespace std;
int main()
{
    double ans = heuman_lambda(sqrt(0.08), 1.5707963267948966);
    cout << ans << endl;
    return 0;
}
