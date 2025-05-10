(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

(* In the Wolfram function repository, m is defined as kc^2. *)

precision = 17;
xRange = DeleteCases[Subdivide[-2, 2, 30], 0];
kcRange = DeleteCases[Subdivide[-2, 2, 30], 0];

(* el1 test data - Complete elliptic integral of the first kind *)
(* Regular data points *)
kcForEl1 = Delete[kcRange, {{1}, {-1}}]; (* Remove extreme values *)
el1 = Flatten[Table[{
  SetPrecision[x, precision],
  SetPrecision[kc, precision], 
  N[ResourceFunction["BulirschEL1"][x, kc*kc], precision]
}, {x, xRange}, {kc, kcForEl1}], 1];
Export["../data/wolfram/el1_data.csv", el1, "CSV"];

(* el1 limit test data - approaching kc->0 *)
kcNearZero = SetPrecision[{0.1, 0.05, 0.01, 0.005, 0.001}, precision];
el1Limit = Flatten[Table[{
  SetPrecision[x, precision],
  kc, 
  N[ResourceFunction["BulirschEL1"][x, kc*kc], precision]
}, {x, xRange}, {kc, kcNearZero}], 1];
Export["../data/wolfram/el1_limit.csv", el1Limit, "CSV"];

(* el2 test data - Complete elliptic integral of the second kind *)
(* Regular data points *)
kcForEl2 = Delete[kcRange, {{1}, {-1}}]; (* Remove extreme values *)
xRange = DeleteCases[Subdivide[-2, 2, 10], 0];
kcRange = DeleteCases[Subdivide[-2, 2, 20], 0];
aRange = Subdivide[-2, 2, 5];
bRange = Subdivide[-2, 2, 5];

el2 = Flatten[Table[{
  SetPrecision[x, precision],
  SetPrecision[kc, precision], 
  SetPrecision[a, precision],
  SetPrecision[b, precision],
  N[ResourceFunction["BulirschEL2"][x, kc * kc, a, b], precision]
}, {x, xRange}, {kc, kcForEl2}, {a, aRange}, {b, bRange}], 3];
el2 = Select[el2, Element[#[[5]], Reals] &];
Export["../data/wolfram/el2_data.csv", el2, "CSV"];

(* el2 limit test data - approaching kc->0 *)
el2Limit = Flatten[Table[{
  SetPrecision[x, precision],
  kc,
  SetPrecision[a, precision],
  SetPrecision[b, precision],
  N[ResourceFunction["BulirschEL2"][x, kc * kc, a, b], precision]
}, {x, xRange}, {kc, kcNearZero}, {a, aRange}, {b, bRange}], 3];
el2Limit = Select[el2Limit, Element[#[[5]], Reals] &];
Export["../data/wolfram/el2_limit.csv", el2Limit, "CSV"];

(* el3 normal test data *)
xRange = DeleteCases[Subdivide[-2, 2, 20], 0];
kcRange = DeleteCases[Subdivide[-2, 2, 20], 0];
pRange = DeleteCases[Subdivide[0.01, 1, 10], 0];
pRangePV = Subdivide[-1, -0.01, 10];

(* Avoid: kc = 0, 1+px^2 = 0, and |kc| > 1 for p < 0 *)
validKcPX = Flatten[Table[
  With[{x = SetPrecision[x, precision], kc = SetPrecision[kc, precision], p = SetPrecision[p, precision]},
    If[kc != 0 && (1 + p*x^2) != 0 && (p >= 0 || Abs[kc] <= 1), 
      {x, kc, p}, 
      Nothing
    ]
  ],
  {x, xRange}, {kc, kcRange}, {p, pRange}
], 2];

el3 = Table[
  With[{x = data[[1]], kc = data[[2]], p = data[[3]]},
    {x, kc, p, N[ResourceFunction["BulirschEL3"][x, kc * kc, p], precision]}
  ],
  {data, validKcPX}
];
el3 = Select[el3, Element[#[[4]], Reals] &];
Export["../data/wolfram/el3_data.csv", el3, "CSV"];

(* el3 with Cauchy principal value when 1 + px^2 < 0 *)
validKcPXPV = Flatten[Table[
  With[{x = SetPrecision[x, precision], kc = SetPrecision[kc, precision], p = SetPrecision[p, precision]},
    If[kc != 0 && (1 + p*x^2) < 0 && (p >= 0 || Abs[kc] <= 1), 
      {x, kc, p}, 
      Nothing
    ]
  ],
  {x, xRange}, {kc, kcRange}, {p, pRangePV}
], 2];

el3PV = Table[
  With[{x = data[[1]], kc = data[[2]], p = data[[3]]},
    {x, kc, p, N[ResourceFunction["BulirschEL3"][x, kc * kc, p], precision]}
  ],
  {data, validKcPXPV}
];
el3PV = Select[el3PV, Element[#[[4]], Reals] &];
Export["../data/wolfram/el3_pv.csv", el3PV, "CSV"];

(* el3 limit test data - approaching critical values *)
(* Near kc = 0 *)
kcNearZero = SetPrecision[{0.1, 0.05, 0.01, 0.005, 0.001}, precision];
xForLimit = SetPrecision[{0.5, 1.0, 1.5}, precision];
pForLimit = SetPrecision[{0.1, 0.5, 1.0}, precision];

el3LimitKc = Flatten[Table[
  {x, kc, p, N[ResourceFunction["BulirschEL3"][x, kc * kc, p], precision]},
  {x, xForLimit}, {kc, kcNearZero}, {p, pForLimit}
], 2];
el3LimitKc = Select[el3LimitKc, Element[#[[4]], Reals] &];
Export["../data/wolfram/el3_limit.csv", el3LimitKc, "CSV"];