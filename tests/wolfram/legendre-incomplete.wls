(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

epsilon = 2.22*^-16;

mu = 1 * epsilon;

PHI = Subdivide[0, Pi / 2, 25];

M = Subdivide[0, 2, 50];

MNeg = Subdivide[-1, 0, 25];

pairs = Flatten[Table[{phi, m}, {phi, PHI}, {m, M}], 1];

normalPairs =
    Select[
        pairs
        ,
        With[{x = 1 - #[[2]] * Sin[#[[1]]] ^ 2},
            x >= mu
        ]&
    ];

(* Incomplete Elliptic Integrals of the First Kind *)

ellipF =
    Table[
        With[{phi = pair[[1]], m = pair[[2]]},
            {phi, m, EllipticF[phi, m]}
        ]
        ,
        {pair, normalPairs}
    ];

Export["../data/wolfram/ellipf_data.csv", N[ellipF], "CSV"];

ellipFNeg = Flatten[Table[{phi, m, EllipticF[phi, m]}, {phi, PHI}, {m,
     MNeg}], 1]

Export["../data/wolfram/ellipf_neg.csv", N[ellipFNeg], "CSV"];

(* Incomplete Elliptic Integrals of the Second Kind *)

ellipEInc =
    Table[
        With[{phi = pair[[1]], m = pair[[2]]},
            {phi, m, EllipticE[phi, m]}
        ]
        ,
        {pair, normalPairs}
    ];

Export["../data/wolfram/ellipeinc_data.csv", N[ellipEInc], "CSV"];

ellipEIncNeg = Flatten[Table[{phi, m, EllipticE[phi, m]}, {phi, PHI},
     {m, MNeg}], 1]

Export["../data/wolfram/ellipeinc_neg.csv", N[ellipEIncNeg], "CSV"];
