(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

(* In the Wolfram function repository, m is defined as kc^2. *)

epsilon = 2.2204460492503131*^-16;

muEL1 = epsilon;

muEL2 = epsilon;

muEL3 = 50 * epsilon;

(* Bulirsch's Incomplete Elliptic integrals of the First Kind *)

xRange = Subdivide[-2, 2, 30];

kcRange = Union[Subdivide[-2, -muEL1, 15], Subdivide[muEL1, 2, 15]];

el1 = Flatten[Table[{x, kc, EllipticF[SetPrecision[ArcTan[x], 30], SetPrecision[
  1 - kc * kc, 30]]}, {x, xRange}, {kc, kcRange}], 1];

Export["../data/wolfram/el1_data.csv", N[el1], "CSV"];

(* Bulirsch's Incomplete Elliptic integrals of the Second Kind *)

xRange = Subdivide[-1, 1, 10];

kcRange = Union[Subdivide[-1, -muEL2, 10], Subdivide[muEL2, 1, 10]];

aRange = Subdivide[-1, 1, 5];

bRange = Subdivide[-1, 1, 5];

el2 = Select[Flatten[Table[{x, kc, a, b, ResourceFunction["BulirschEL2"
  ][SetPrecision[x, 30], SetPrecision[kc * kc, 30], SetPrecision[a, 30],
   SetPrecision[b, 30]]}, {x, xRange}, {kc, kcRange}, {a, aRange}, {b, 
  bRange}], 3], Im[#[[5]]] == 0&];

Export["../data/wolfram/el2_data.csv", N[el2], "CSV"];

(* Bulirsch's Incomplete Elliptic integrals of the Third Kind - Normal domain *)

xRange = Union[Subdivide[-1 + muEL3, 1 - muEL3, 10], {0}];

kcRange = Union[Subdivide[-1.5, -muEL3, 10], Subdivide[muEL3, 1.5, 10
  ]];

pRange = Subdivide[muEL3, 1.0, 10];

pRangePV = Subdivide[-3.0, -1.1, 10];

EL3[x_, kc_, p_] :=
  EllipticPi[SetPrecision[1 - p, 30], SetPrecision[ArcTan[x], 30], SetPrecision[
    1 - kc * kc, 30]];

validKcPX =
  Select[
    Flatten[Table[{x, kc, p}, {x, xRange}, {kc, kcRange}, {p, pRange}
      ], 2]
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

Export["../data/wolfram/el3_data.csv", N[el3], "CSV"];

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

Export["../data/wolfram/el3_pv.csv", N[el3PV], "CSV"];
