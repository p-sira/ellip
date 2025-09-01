(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

(* ========== f32 ========== *)

epsilon = 1.19*^-7;

mu = 1 * epsilon;

M = Subdivide[0, 1 - mu, 1000];

MNeg = Subdivide[-5, -mu, 1000];

(* ----- Complete Elliptic Integrals -----*)

(* Complete Elliptic Integrals of the First Kind *)

ellipK = Table[{m, EllipticK[SetPrecision[m, 30]]}, {m, M}];

Export["../data/f32/wolfram/ellipk_data.csv", N[ellipK], "CSV"];

ellipKNeg = Table[{m, EllipticK[SetPrecision[m, 30]]}, {m, MNeg}];

Export["../data/f32/wolfram/ellipk_neg.csv", N[ellipKNeg], "CSV"];

(* Complete Elliptic Integrals of the Second Kind *)

ellipE = Table[{m, EllipticE[SetPrecision[m, 30]]}, {m, M}];

Export["../data/f32/wolfram/ellipe_data.csv", N[ellipE], "CSV"];

ellipENeg = Table[{m, EllipticE[SetPrecision[m, 30]]}, {m, MNeg}];

Export["../data/f32/wolfram/ellipe_neg.csv", N[ellipENeg], "CSV"];

(* Complete Elliptic Integrals of Legendre's Type *)

M = Subdivide[mu, 1 - mu, 1000];

EllipticD[m_] :=
    (EllipticK[m] - EllipticE[m]) / m;

ellipD = Table[{m, EllipticD[SetPrecision[m, 30]]}, {m, M}];

Export["../data/f32/wolfram/ellipd_data.csv", N[ellipD], "CSV"];

ellipDNeg = Table[{m, EllipticD[SetPrecision[m, 30]]}, {m, MNeg}];

Export["../data/f32/wolfram/ellipd_neg.csv", N[ellipDNeg], "CSV"];

muPV = 50 * epsilon;

mNormal = Subdivide[0, 1 - mu, 50];

ns = Subdivide[0, 1 - mu, 20];

(* Elliptic integrals of the Third Kind - Normal domain *)

ellipPiNormal = Flatten[Table[{n, m, EllipticPi[SetPrecision[n, 30], 
    SetPrecision[m, 30]]}, {n, ns}, {m, mNormal}], 1];

ellipPiNormal = Select[ellipPiNormal, Element[#[[3]], Reals]&];

Export["../data/f32/wolfram/ellippi_data.csv", N[ellipPiNormal], "CSV"
    ];

(* Elliptic integrals of the Third Kind - Negative m domain *)

mNeg = Subdivide[-5, -0.001, 50];

ellipPiNeg = Flatten[Table[{n, m, EllipticPi[SetPrecision[n, 30], SetPrecision[
    m, 30]]}, {n, ns}, {m, mNeg}], 1];

ellipPiNeg = Select[ellipPiNeg, Element[#[[3]], Reals]&];

Export["../data/f32/wolfram/ellippi_neg.csv", N[ellipPiNeg], "CSV"];

(* Elliptic integrals of the Third Kind - PV domain *)

nPV = Subdivide[1 + muPV, 10, 20];

mPV = Subdivide[0, 1 - muPV, 25];

ellipPiPV = Flatten[Table[{n, m, EllipticPi[SetPrecision[n, 30], SetPrecision[
    m, 30]]}, {n, nPV}, {m, mPV}], 1];

ellipPiPV = Select[ellipPiPV, Element[#[[3]], Reals]&];

Export["../data/f32/wolfram/ellippi_pv.csv", N[ellipPiPV], "CSV"];

(* ----- Incomplete Elliptic Integrals ----- *)

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
            {phi, m, EllipticF[SetPrecision[phi, 30], SetPrecision[m,
                 30]]}
        ]
        ,
        {pair, normalPairs}
    ];

Export["../data/f32/wolfram/ellipf_data.csv", N[ellipF], "CSV"];

ellipFNeg = Flatten[Table[{phi, m, EllipticF[SetPrecision[phi, 30], SetPrecision[
    m, 30]]}, {phi, PHI}, {m, MNeg}], 1]

Export["../data/f32/wolfram/ellipf_neg.csv", N[ellipFNeg], "CSV"];

(* Incomplete Elliptic Integrals of the Second Kind *)

ellipEInc =
    Table[
        With[{phi = pair[[1]], m = pair[[2]]},
            {phi, m, EllipticE[SetPrecision[phi, 30], SetPrecision[m,
                 30]]}
        ]
        ,
        {pair, normalPairs}
    ];

Export["../data/f32/wolfram/ellipeinc_data.csv", N[ellipEInc], "CSV"];

ellipEIncNeg = Flatten[Table[{phi, m, EllipticE[SetPrecision[phi, 30],
     SetPrecision[m, 30]]}, {phi, PHI}, {m, MNeg}], 1]

Export["../data/f32/wolfram/ellipeinc_neg.csv", N[ellipEIncNeg], "CSV"
    ];

(* Incomplete Elliptic Integrals of Legendre's Type *)

EllipticD[phi_, m_] :=
    Module[{c},
        Which[
            phi == 0,
                0
            ,
            True,
                c = Csc[phi] ^ 2;
                CarlsonRD[c - 1, c - m, c] / 3
        ]
    ];

ellipDInc =
    Table[
        With[{phi = pair[[1]], m = pair[[2]]},
            {phi, m, EllipticD[SetPrecision[phi, 30], SetPrecision[m,
                 30]]}
        ]
        ,
        {pair, normalPairs}
    ];

Export["../data/f32/wolfram/ellipdinc_data.csv", N[ellipDInc], "CSV"];

ellipDIncNeg = Flatten[Table[{phi, m, EllipticD[SetPrecision[phi, 30],
     SetPrecision[m, 30]]}, {phi, PHI}, {m, MNeg}], 1]

Export["../data/f32/wolfram/ellipdinc_neg.csv", N[ellipDIncNeg], "CSV"
    ];

(* Incomplete Elliptic Integrals of the Third Kind *)

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
                {phi, n, m, EllipticPi[SetPrecision[n, 30], SetPrecision[
                    phi, 30], SetPrecision[m, 30]]}
            ]
            ,
            {triplet, normalTriplets}
        ]
        ,
        Im[#[[4]]] == 0&
    ];

Export["../data/f32/wolfram/ellippiinc_data.csv", N[ellipPiInc], "CSV"
    ];

(* Elliptic integrals of the Third Kind - Negative m domain *)

ellipPiIncNeg =
    Flatten[
        Table[
            With[{phi = phiVal, m = mVal, n = nVal},
                {phi, n, m, EllipticPi[SetPrecision[n, 30], SetPrecision[
                    phi, 30], SetPrecision[m, 30]]}
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

Export["../data/f32/wolfram/ellippiinc_neg.csv", N[ellipPiIncNeg], "CSV"
    ];

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
                {phi, n, m, EllipticPi[SetPrecision[n, 30], SetPrecision[
                    phi, 30], SetPrecision[m, 30]]}
            ]
            ,
            {triplet, pvTriplets}
        ]
        ,
        Im[#[[4]]] == 0&
    ];

Export["../data/f32/wolfram/ellippiinc_pv.csv", N[ellipPiIncPV], "CSV"
    ];
