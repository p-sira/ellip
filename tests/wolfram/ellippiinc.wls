(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

epsilon = 2.22*^-16;

mu = 50 * epsilon;

PHI = Subdivide[0, Pi / 2, 5];

mNormal = Subdivide[0, 2, 35];

mNeg = Subdivide[-1, 0, 15];

mPV = Subdivide[0, 1, 80];

NVals = Subdivide[-1, 2, 5];

(* Elliptic integrals of the Third Kind - Normal domain *)

normalTriplets =
    Select[
        Flatten[Table[{phi, m, n}, {phi, PHI}, {m, mNormal}, {n, NVals
            }], 2]
        ,
        With[{phi = #[[1]], m = #[[2]], n = #[[3]], s2 = Sin[#[[1]]] 
            ^ 2},
            m * s2 < 1 - mu || n * s2 != 1
        ]&
    ];

ellipPiInc =
    Select[
        Table[
            With[{phi = triplet[[1]], m = triplet[[2]], n = triplet[[
                3]]},
                {phi, n, m, EllipticPi[n, phi, m]}
            ]
            ,
            {triplet, normalTriplets}
        ]
        ,
        Im[#[[4]]] == 0&
    ];

Export["../data/wolfram/ellippiinc_data.csv", N[ellipPiInc], "CSV"];

(* Elliptic integrals of the Third Kind - Negative m domain *)

ellipPiIncNeg =
    Flatten[
        Table[
            With[{phi = phiVal, m = mVal, n = nVal},
                {phi, n, m, EllipticPi[n, phi, m]}
            ]
            ,
            {phiVal, PHI}
            ,
            {mVal, mNeg}
            ,
            {nVal, NVals}
        ]
        ,
        2
    ];

Export["../data/wolfram/ellippiinc_neg.csv", N[ellipPiIncNeg], "CSV"];

(* Elliptic integrals of the Third Kind - PV domain *)

pvTriplets =
    Select[
        Flatten[Table[{phi, m, n}, {phi, PHI}, {m, mPV}, {n, NVals}],
             2]
        ,
        With[{phi = #[[1]], m = #[[2]], n = #[[3]], s2 = Sin[#[[1]]] 
            ^ 2},
            n * s2 > 1 + mu && m * s2 < 1 - mu
        ]&
    ];

ellipPiIncPV =
    Select[
        Table[
            With[{phi = triplet[[1]], m = triplet[[2]], n = triplet[[
                3]]},
                {phi, n, m, EllipticPi[n, phi, m]}
            ]
            ,
            {triplet, pvTriplets}
        ]
        ,
        Im[#[[4]]] == 0&
    ];

Export["../data/wolfram/ellippiinc_pv.csv", N[ellipPiIncPV], "CSV"];
