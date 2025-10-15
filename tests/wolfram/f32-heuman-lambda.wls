(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

epsilon = 1.1920929*^-7;

mu = 1 * epsilon;

PHI = Subdivide[0, 2 * Pi, 40];

M = Subdivide[mu, 1 - mu, 25];

HeumanLambda[phi_, m_] :=
    Module[{mc},
        mc = 1 - m;
        EllipticF[phi, mc] / EllipticK[mc] + (2 / Pi) * EllipticK[m] 
            * JacobiZeta[phi, mc]
    ]

lambda = Flatten[Table[{phi, m, HeumanLambda[SetPrecision[phi, 30], SetPrecision[
    m, 30]]}, {phi, PHI}, {m, M}], 1];

Export["../data/f32/wolfram/heuman_lambda_data.csv", N[lambda], "CSV"
    ];
