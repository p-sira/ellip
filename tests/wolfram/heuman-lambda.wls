(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

epsilon = 2.22*^-16;

mu = 1 * epsilon;

PHI = Subdivide[0, 2 * Pi, 40];

M = Subdivide[mu, 1 - mu, 25];

HeumanLambda[phi_, m_] :=
    Module[{smc},
        smc = Sqrt[1 - m];
        EllipticF[phi, smc] / EllipticK[smc] + (2 / Pi) * EllipticK[m
            ] * JacobiZeta[phi, smc]
    ]

lambda = Flatten[Table[{phi, m, HeumanLambda[SetPrecision[phi, 30], SetPrecision[
    m, 30]]}, {phi, PHI}, {m, M}], 1];

Export["../data/wolfram/heuman_lambda_data.csv", N[lambda], "CSV"];
