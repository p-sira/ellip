(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

(* In the Wolfram function repository, m is defined as kc^2. *)

epsilon = 2.22*^-16;

mu = 1 * epsilon;

kcRange = Union[Subdivide[-2, -mu, 5], Subdivide[mu, 2, 5]];

pRange = Subdivide[mu, 1, 5];

pRangePV = Subdivide[-1, -mu, 5];

aRange = Subdivide[-2, 2, 4];

bRange = Subdivide[-2, 2, 4];

cel = Select[Flatten[Table[{kc, p, a, b, ResourceFunction["BulirschCEL"
    ][kc * kc, p, a, b]}, {kc, kcRange}, {p, pRange}, {a, aRange}, {b, bRange
    }], 3], Element[#[[5]], Reals]&];

celPV = Select[Flatten[Table[{kc, p, a, b, ResourceFunction["BulirschCEL"
    ][kc * kc, p, a, b]}, {kc, kcRange}, {p, pRangePV}, {a, aRange}, {b, 
    bRange}], 3], Element[#[[5]], Reals]&];

Export["../data/wolfram/cel_data.csv", N[cel], "CSV"];

Export["../data/wolfram/cel_pv.csv", N[celPV], "CSV"];
