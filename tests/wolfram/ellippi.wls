(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

precision = 15;

mNormal = Subdivide[0, 0.999, 100];

ns = Subdivide[0, 0.9, 50];

ellipPiNormal = Flatten[Table[{SetPrecision[n, precision], SetPrecision[
    m, precision], N[EllipticPi[n, m], precision]}, {n, ns}, {m, mNormal}
    ], 1];

Export["../data/wolfram/ellippi_data.csv", ellipPiNormal, "CSV"];

mNeg = Subdivide[-5, -0.001, 100];

ellipPiNeg = Flatten[Table[{SetPrecision[n, precision], SetPrecision[
    m, precision], N[EllipticPi[n, m], precision]}, {n, ns}, {m, mNeg}], 
    1];

Export["../data/wolfram/ellippi_neg.csv", ellipPiNeg, "CSV"];

nPV = Subdivide[1.001, 10, 50];

mPV = Subdivide[0, 0.99, 50];

ellipPiPV = Flatten[Table[{SetPrecision[n, precision], SetPrecision[m,
     precision], N[EllipticPi[n, m], precision]}, {n, nPV}, {m, mPV}], 1];

Export["../data/wolfram/ellippi_pv.csv", ellipPiPV, "CSV"];

mLimit = Subdivide[0.999, 1, 50];

ellipPiLimit = Flatten[Table[{SetPrecision[n, precision], SetPrecision[
    m, precision], N[EllipticPi[n, m], precision]}, {n, ns}, {m, mLimit}],
     1];

Export["../data/wolfram/ellippi_limit.csv", ellipPiLimit, "CSV"];
