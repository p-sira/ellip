(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

precision = 17;

epsilon = 1*^-3;

(* Singularity condition: 1 - m sin^2(phi) *)

MLim = Subdivide[0.1, 2, 10];

(* Sample phi near the limit given m *)

limitPairs = Reap[
        Do[
            If[m > epsilon,
                Module[{s2Low, s2High, phiLow, phiHigh, phiSample},
                    s2Low = (1 - epsilon) / m;
                    s2High = (1 + epsilon) / m;
                    If[0 < s2Low < 1 && 0 < s2High < 1,
                        phiLow = ArcSin[Sqrt[s2Low]];
                        phiHigh = ArcSin[Sqrt[s2High]];
                        phiSample = Subdivide[phiLow, phiHigh, 10];
                        Do[Sow[{phi, m}], {phi, phiSample}];
                    ];
                ]
            ]
            ,
            {m, MLim}
        ]
    ][[2, 1]];

(* Filter in only params that yield answers within the real domain *)

limitPairs =
    Select[
        limitPairs
        ,
        With[{phi = #[[1]], m = #[[2]]},
            0 < 1 - m Sin[phi] ^ 2 < epsilon
        ]&
    ];

PHI = Subdivide[0, Pi / 2, 50];

M = Subdivide[0, 2, 100];

MNeg = Subdivide[-1, 0, 50];

pairs = Flatten[Table[{phi, m}, {phi, PHI}, {m, M}], 1];

normalPairs =
    Select[
        pairs
        ,
        With[{x = 1 - #[[2]] * Sin[#[[1]]] ^ 2},
            x >= epsilon
        ]&
    ];

ellipF =
    Table[
        With[{phi = pair[[1]], m = pair[[2]]},
            {SetPrecision[phi, precision], SetPrecision[m, precision],
                 N[EllipticF[phi, m], precision]}
        ]
        ,
        {pair, normalPairs}
    ];

Export["../data/wolfram/ellipf_data.csv", ellipF, "CSV"];

ellipFNeg = Flatten[Table[{SetPrecision[phi, precision], SetPrecision[
    m, precision], N[EllipticF[phi, m], precision]}, {phi, PHI}, {m, MNeg
    }], 1]

Export["../data/wolfram/ellipf_neg.csv", ellipFNeg, "CSV"];

ellipFLimit =
    Select[
        Table[
            With[{phi = pair[[1]], m = pair[[2]], result = N[EllipticF[
                pair[[1]], pair[[2]]], precision]},
                {SetPrecision[phi, precision], SetPrecision[m, precision
                    ], result}
            ]
            ,
            {pair, limitPairs}
        ]
        ,
        (* Keep only cases with real answers *)
        Im[#[[3]]] == 0&
    ];

Export["../data/wolfram/ellipf_limit.csv", ellipFLimit, "CSV"];

ellipEInc =
    Table[
        With[{phi = pair[[1]], m = pair[[2]]},
            {SetPrecision[phi, precision], SetPrecision[m, precision],
                 N[EllipticE[phi, m], precision]}
        ]
        ,
        {pair, normalPairs}
    ];

Export["../data/wolfram/ellipeinc_data.csv", ellipEInc, "CSV"];

ellipEIncNeg = Flatten[Table[{SetPrecision[phi, precision], SetPrecision[
    m, precision], N[EllipticE[phi, m], precision]}, {phi, PHI}, {m, MNeg
    }], 1]

Export["../data/wolfram/ellipeinc_neg.csv", ellipEIncNeg, "CSV"];

ellipEIncLimit =
    Select[
        Table[
            With[{phi = pair[[1]], m = pair[[2]], result = N[EllipticE[
                pair[[1]], pair[[2]]], precision]},
                {SetPrecision[phi, precision], SetPrecision[m, precision
                    ], result}
            ]
            ,
            {pair, limitPairs}
        ]
        ,
        (* Keep only cases with real answers *)
        Im[#[[3]]] == 0&
    ];

Export["../data/wolfram/ellipeinc_limit.csv", ellipEIncLimit, "CSV"];
