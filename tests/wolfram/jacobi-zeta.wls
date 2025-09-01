(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

epsilon = 2.22*^-16;

mu = 1 * epsilon;

PHI = Subdivide[0, 2 * Pi, 40];

M = Subdivide[0, 1, 25];

zeta = Flatten[Table[{phi, m, JacobiZeta[SetPrecision[phi, 30], SetPrecision[
    m, 30]]}, {phi, PHI}, {m, M}], 1];

Export["../data/wolfram/jacobi_zeta_data.csv", N[zeta], "CSV"];

MNeg = Subdivide[-2, 0, 25];

zetaNeg = Flatten[Table[{phi, m, JacobiZeta[SetPrecision[phi, 30], SetPrecision[
    m, 30]]}, {phi, PHI}, {m, MNeg}], 1];

Export["../data/wolfram/jacobi_zeta_neg.csv", N[zetaNeg], "CSV"];
