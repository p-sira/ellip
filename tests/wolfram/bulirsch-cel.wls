(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

(* In the Wolfram function repository, m is defined as kc^2. *)

precision = 17;

kcRange = DeleteCases[Subdivide[-2, 2, 20], 0];

pRange = DeleteCases[Subdivide[0.01, 1, 10], 0];

pRangePV = Subdivide[-1, -0.01, 10];

aRange = Subdivide[-2, 2, 5];

bRange = Subdivide[-2, 2, 5];

cel = Select[Flatten[Table[{SetPrecision[kc, precision], SetPrecision[
    p, precision], SetPrecision[a, precision], SetPrecision[b, precision],
     N[ResourceFunction["BulirschCEL"][kc * kc, p, a, b], precision]}, {kc,
     kcRange}, {p, pRange}, {a, aRange}, {b, bRange}], 3], Element[#[[5]],
     Reals]&];

celPV = Select[Flatten[Table[{SetPrecision[kc, precision], SetPrecision[
    p, precision], SetPrecision[a, precision], SetPrecision[b, precision],
     N[ResourceFunction["BulirschCEL"][kc * kc, p, a, b], precision]}, {kc,
     kcRange}, {p, pRangePV}, {a, aRange}, {b, bRange}], 3], Element[#[[5
    ]], Reals]&];

Export["../data/wolfram/cel_data.csv", cel, "CSV"];

Export["../data/wolfram/cel_pv.csv", celPV, "CSV"];
