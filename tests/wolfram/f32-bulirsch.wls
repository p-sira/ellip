(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

(* ========== f32 ========== *)

(* In the Wolfram function repository, m is defined as kc^2. *)

epsilon = 1.19*^-7;

(* ----- Bulirsch's Complete Elliptic integrals ----- *)

muCel = 1 * epsilon;

muCel1 = 1 * epsilon;

muCel2 = 1 * epsilon;

(* Bulirsch's General Complete Elliptic integrals *)

kcRange = Union[Subdivide[-2, -muCel, 5], Subdivide[muCel, 2, 5]];

pRange = Subdivide[muCel, 1, 5];

pRangePV = Subdivide[-1, -muCel, 5];

aRange = Subdivide[-2, 2, 4];

bRange = Subdivide[-2, 2, 4];

cel = Select[Flatten[Table[{kc, p, a, b, ResourceFunction["BulirschCEL"
    ][SetPrecision[kc * kc, 30], SetPrecision[p, 30], SetPrecision[a, 30],
     SetPrecision[b, 30]]}, {kc, kcRange}, {p, pRange}, {a, aRange}, {b, 
    bRange}], 3], Element[#[[5]], Reals]&];

celPV = Select[Flatten[Table[{kc, p, a, b, ResourceFunction["BulirschCEL"
    ][SetPrecision[kc * kc, 30], SetPrecision[p, 30], SetPrecision[a, 30],
     SetPrecision[b, 30]]}, {kc, kcRange}, {p, pRangePV}, {a, aRange}, {b,
     bRange}], 3], Element[#[[5]], Reals]&];

Export["../data/f32/wolfram/cel_data.csv", N[cel], "CSV"];

Export["../data/f32/wolfram/cel_pv.csv", N[celPV], "CSV"];

(* Bulirsch's Complete Elliptic integrals of the First Kind *)

kcRange = Union[Subdivide[-2, -muCel1, 50], Subdivide[muCel1, 2, 50]];

cel1 = Select[Table[{kc, EllipticK[SetPrecision[1 - kc * kc, 30]]}, {
    kc, kcRange}], Element[#[[2]], Reals]&];

Export["../data/f32/wolfram/cel1_data.csv", N[cel1], "CSV"];

(* Bulirsch's Complete Elliptic integrals of the Second Kind *)

kcRange = Union[Subdivide[-2, -muCel2, 10], Subdivide[muCel2, 2, 10]];

cel2 = Select[Flatten[Table[{kc, a, b, ResourceFunction["BulirschCEL"
    ][SetPrecision[kc * kc, 30], 1, SetPrecision[a, 30], SetPrecision[b, 
    30]]}, {kc, kcRange}, {a, aRange}, {b, bRange}], 2], Element[#[[4]], 
    Reals]&];

Export["../data/f32/wolfram/cel2_data.csv", N[cel2], "CSV"];

(* ----- Bulirsch's Incomplete Elliptic integrals ----- *)

muEL1 = epsilon;

muEL2 = epsilon;

muEL3 = 50 * epsilon;

(* Bulirsch's Incomplete Elliptic integrals of the First Kind *)

xRange = Subdivide[-2, 2, 30];

kcRange = Union[Subdivide[-2, -muEL1, 15], Subdivide[muEL1, 2, 15]];

el1 = Flatten[Table[{x, kc, EllipticF[SetPrecision[ArcTan[x], 30], SetPrecision[
    1 - kc * kc, 30]]}, {x, xRange}, {kc, kcRange}], 1];

Export["../data/f32/wolfram/el1_data.csv", N[el1], "CSV"];

(* Bulirsch's Incomplete Elliptic integrals of the Second Kind *)

xRange = Subdivide[-1, 1, 10];

kcRange = Union[Subdivide[-1, -muEL2, 10], Subdivide[muEL2, 1, 10]];

aRange = Subdivide[-1, 1, 5];

bRange = Subdivide[-1, 1, 5];

el2 = Select[Flatten[Table[{x, kc, a, b, ResourceFunction["BulirschEL2"
    ][SetPrecision[x, 30], SetPrecision[kc * kc, 30], SetPrecision[a, 30],
     SetPrecision[b, 30]]}, {x, xRange}, {kc, kcRange}, {a, aRange}, {b, 
    bRange}], 3], Im[#[[5]]] == 0&];

Export["../data/f32/wolfram/el2_data.csv", N[el2], "CSV"];

(* Bulirsch's Incomplete Elliptic integrals of the Third Kind - Normal domain *)

xRange = Union[Subdivide[-1 + muEL3, 1 - muEL3, 10], {0}];

kcRange = Union[Subdivide[-1.5, -muEL3, 10], Subdivide[muEL3, 1.5, 10
    ]];

pRange = Subdivide[muEL3, 1.0, 10];

pRangePV = Subdivide[-3.0, -1.1, 10];

EL3[x_, kc_, p_] :=
    EllipticPi[SetPrecision[1 - p, 30], SetPrecision[ArcTan[x], 30], 
        SetPrecision[1 - kc * kc, 30]];

validKcPX =
    Select[
        Flatten[Table[{x, kc, p}, {x, xRange}, {kc, kcRange}, {p, pRange
            }], 2]
        ,
        With[{x = #[[1]], kc = #[[2]], p = #[[3]]},
            kc != 0 && (1 + p * x^2) > muEL3
        ]&
    ];

el3 =
    Select[
        Table[
            With[{x = data[[1]], kc = data[[2]], p = data[[3]]},
                {x, kc, p, EL3[x, kc, p]}
            ]
            ,
            {data, validKcPX}
        ]
        ,
        Im[#[[4]]] == 0&
    ];

Export["../data/f32/wolfram/el3_data.csv", N[el3], "CSV"];

(* Bulirsch's Incomplete Elliptic integrals of the Third Kind - PV domain (1 + px^2 < 0) *)

validKcPXPV =
    Select[
        Flatten[Table[{x, kc, p}, {x, xRange}, {kc, kcRange}, {p, pRangePV
            }], 2]
        ,
        With[{x = #[[1]], kc = #[[2]], p = #[[3]]},
            kc != 0 && (1 + p * x^2) < -muEL3
        ]&
    ];

el3PV =
    Select[
        Table[
            With[{x = data[[1]], kc = data[[2]], p = data[[3]]},
                {x, kc, p, EL3[x, kc, p]}
            ]
            ,
            {data, validKcPXPV}
        ]
        ,
        Im[#[[4]]] == 0&
    ];

Export["../data/f32/wolfram/el3_pv.csv", N[el3PV], "CSV"];
