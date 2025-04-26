(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

precision = 15;

epsilon = 1*^-4;

(* Sample domain limits *)

PHI = Subdivide[0, Pi / 2, 20];

M = Subdivide[0, 2, 40];

MNeg = Subdivide[-1, 0, 40];

NVals = Subdivide[-1, 2, 10];

normalTriplets =
    Select[
        Flatten[Table[{phi, m, n}, {phi, PHI}, {m, M}, {n, NVals}], 2
            ]
        ,
        With[{phi = #[[1]], m = #[[2]], n = #[[3]], s2 = Sin[#[[1]]] 
            ^ 2},
            m * s2 < 1 || n * s2 != 1
        ]&
    ];

MLim = Subdivide[1.001, 3, 10]

epsilon = 1*^-4;

limitTriplets = Reap[
        Do[
            If[m > epsilon,
                Module[{s2Low, s2High, phiLow, phiHigh, phiSample},
                    s2Low = (1 - epsilon) / m;
                    s2High = (1 + epsilon) / m;
                    If[0 < s2Low < 1 && 0 < s2High < 1,
                        phiLow = ArcSin[Sqrt[s2Low]];
                        phiHigh = ArcSin[Sqrt[s2High]];
                        phiSample = Subdivide[phiLow, phiHigh, 5];
                        Do[Sow[{phi, m, n}], {phi, phiSample}, {n, NVals
                            }];
                    ];
                ]
            ]
            ,
            {m, MLim}
        ]
    ][[2, 1]];

(* Normal domain *)

ellipPiInc =
    Select[
        Table[
            With[{phi = triplet[[1]], m = triplet[[2]], n = triplet[[
                3]]},
                {SetPrecision[phi, precision], SetPrecision[n, precision
                    ], SetPrecision[m, precision], N[EllipticPi[n, phi, m], precision]}
            ]
            ,
            {triplet, normalTriplets}
        ]
        ,
        Im[#[[4]]] == 0&
    ];

(* Negative m values *)

ellipPiIncNeg =
    Flatten[
        Table[
            With[{phi = phiVal, m = mVal, n = nVal},
                {SetPrecision[phi, precision], SetPrecision[n, precision
                    ], SetPrecision[m, precision], N[EllipticPi[n, phi, m], precision]}
            ]
            ,
            {phiVal, PHI}
            ,
            {mVal, MNeg}
            ,
            {nVal, NVals}
        ]
        ,
        2
    ];

(* Limit approaching m s2p = 1- *)

ellipPiIncLimit =
    Select[
        Table[
            With[{phi = N[triple[[1]], precision], m = N[triple[[2]],
                 precision], n = N[triple[[3]], precision]},
                {phi, n, m, N[EllipticPi[n, phi, m], precision]}
            ]
            ,
            {triple, limitTriplets}
        ]
        ,
        Im[#[[4]]] == 0&
    ];

Export["../data/wolfram/ellippiinc_data.csv", ellipPiInc, "CSV"];

Export["../data/wolfram/ellippiinc_neg.csv", ellipPiIncNeg, "CSV"];

Export["../data/wolfram/ellippiinc_limit.csv", ellipPiIncLimit, "CSV"
    ];
