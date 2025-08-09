(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

(* Test coverage for all branches (m in [0,1]) *)

idx = Range[0, 19, 1];

M = Union[idx / 20, {1 - 1*^-6}];

(* Complete Elliptic Integrals of the First Kind *)

ellipK = Table[{m, EllipticK[SetPrecision[m, 30]]}, {m, M}];

Export["../data/wolfram/ellipk_cov.csv", N[ellipK], "CSV"];

(* Complete Elliptic Integrals of the Second Kind *)

(* Test coverage for m in [0,1] *)

ellipE = Table[{m, EllipticE[SetPrecision[m, 30]]}, {m, M}];

Export["../data/wolfram/ellipe_cov.csv", N[ellipE], "CSV"];
