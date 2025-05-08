(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

precision = 17;

M = Subdivide[0, 0.999, 1000];

MLimit = Subdivide[0.999, 1, 100];

MNeg = Subdivide[-5, -0.001, 1000];

ellipK = Table[{SetPrecision[m, precision], N[EllipticK[m], precision
    ]}, {m, M}];

Export["../data/wolfram/ellipk_data.csv", ellipK, "CSV"];

ellipKLimit = Table[{SetPrecision[m, precision], N[EllipticK[m], precision
    ]}, {m, MLimit}];

Export["../data/wolfram/ellipk_limit.csv", ellipKLimit, "CSV"];

ellipKNeg = Table[{SetPrecision[m, precision], N[EllipticK[m], precision
    ]}, {m, MNeg}];

Export["../data/wolfram/ellipk_neg.csv", ellipKNeg, "CSV"];

ellipE = Table[{SetPrecision[m, precision], N[EllipticE[m], precision
    ]}, {m, M}];

Export["../data/wolfram/ellipe_data.csv", ellipE, "CSV"];

ellipELimit = Table[{SetPrecision[m, precision], N[EllipticE[m], precision
    ]}, {m, MLimit}];

Export["../data/wolfram/ellipe_limit.csv", ellipELimit, "CSV"];

ellipENeg = Table[{SetPrecision[m, precision], N[EllipticE[m], precision
    ]}, {m, MNeg}];

Export["../data/wolfram/ellipe_neg.csv", ellipENeg, "CSV"];
