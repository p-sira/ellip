(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

(* In the Wolfram function repository, m is defined as kc^2. *)

epsilon = 2.22*^-16;

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
    ][kc * kc, p, a, b]}, {kc, kcRange}, {p, pRange}, {a, aRange}, {b, bRange
    }], 3], Element[#[[5]], Reals]&];

celPV = Select[Flatten[Table[{kc, p, a, b, ResourceFunction["BulirschCEL"
    ][kc * kc, p, a, b]}, {kc, kcRange}, {p, pRangePV}, {a, aRange}, {b, 
    bRange}], 3], Element[#[[5]], Reals]&];

Export["../data/wolfram/cel_data.csv", N[cel], "CSV"];

Export["../data/wolfram/cel_pv.csv", N[celPV], "CSV"];

(* Bulirsch's Complete Elliptic integrals of the First Kind *)

kcRange = Union[Subdivide[-2, -muCel1, 50], Subdivide[muCel1, 2, 50]];

cel1 = Select[Table[{kc, EllipticK[1 - kc * kc]}, {kc, kcRange}], Element[
    #[[2]], Reals]&];

Export["../data/wolfram/cel1_data.csv", N[cel1], "CSV"];

(* Bulirsch's Complete Elliptic integrals of the Second Kind *)

kcRange = Union[Subdivide[-2, -muCel2, 10], Subdivide[muCel2, 2, 10]];

cel2 = Select[Flatten[Table[{kc, a, b, ResourceFunction["BulirschCEL"
    ][kc * kc, 1, a, b]}, {kc, kcRange}, {a, aRange}, {b, bRange}], 2], Element[
    #[[4]], Reals]&];

Export["../data/wolfram/cel2_data.csv", N[cel2], "CSV"];
