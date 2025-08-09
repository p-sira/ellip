(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

epsilon = 2.22*^-16;

(* Mu helps prevent answers approaching infinity, making the comparison fairer. *)

mu = 50 * epsilon

mNormal = Subdivide[0, 1 - mu, 50];

ns = Subdivide[0, 1 - mu, 20];

(* Elliptic integrals of the Third Kind - Normal domain *)

ellipPiNormal = Flatten[Table[{n, m, EllipticPi[SetPrecision[n, 30], 
    SetPrecision[m, 30]]}, {n, ns}, {m, mNormal}], 1];

ellipPiNormal = Select[ellipPiNormal, Element[#[[3]], Reals]&];

Export["../data/wolfram/ellippi_data.csv", N[ellipPiNormal], "CSV"];

(* Elliptic integrals of the Third Kind - Negative m domain *)

mNeg = Subdivide[-5, -0.001, 50];

ellipPiNeg = Flatten[Table[{n, m, EllipticPi[SetPrecision[n, 30], SetPrecision[
    m, 30]]}, {n, ns}, {m, mNeg}], 1];

ellipPiNeg = Select[ellipPiNeg, Element[#[[3]], Reals]&];

Export["../data/wolfram/ellippi_neg.csv", N[ellipPiNeg], "CSV"];

(* Elliptic integrals of the Third Kind - PV domain *)

nPV = Subdivide[1 + mu, 10, 20];

mPV = Subdivide[0, 1 - mu, 25];

ellipPiPV = Flatten[Table[{n, m, EllipticPi[SetPrecision[n, 30], SetPrecision[
    m, 30]]}, {n, nPV}, {m, mPV}], 1];

ellipPiPV = Select[ellipPiPV, Element[#[[3]], Reals]&];

Export["../data/wolfram/ellippi_pv.csv", N[ellipPiPV], "CSV"];
