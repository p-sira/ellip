(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

epsilon = 2.22*^-16;

mu = 1 * epsilon;

M = Subdivide[0, 1 - mu, 1000];

MNeg = Subdivide[-5, -mu, 1000];

(* Complete Elliptic Integrals of the First Kind *)

ellipK = Table[{m, EllipticK[m]}, {m, M}];

Export["../data/wolfram/ellipk_data.csv", ellipK, "CSV"];

ellipKNeg = Table[{m, EllipticK[m]}, {m, MNeg}];

Export["../data/wolfram/ellipk_neg.csv", ellipKNeg, "CSV"];

(* Complete Elliptic Integrals of the Second Kind *)

ellipE = Table[{m, EllipticE[m]}, {m, M}];

Export["../data/wolfram/ellipe_data.csv", ellipE, "CSV"];

ellipENeg = Table[{m, EllipticE[m]}, {m, MNeg}];

Export["../data/wolfram/ellipe_neg.csv", ellipENeg, "CSV"];
