(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

(* ========== f32 ========== *)

epsilon = 1.19*^-7;

muRC = epsilon;

muRD = 50 * epsilon;

muRG = epsilon;

muRF = epsilon;

muRJ = 50 * epsilon;

muRJPV = 1*^10 * epsilon;

muRJPVSmall = 50 * epsilon;

(* Carlson RC - Degenerate symmetric elliptic integral of RF *)

(* Normal domain: x >= 0, y > 0 *)

xRangeRC = Union[Subdivide[muRC, 2, 30], {0}];

yRangeRC = Subdivide[muRC, 2, 30];

rcNormal = Flatten[Table[{x, y, CarlsonRC[SetPrecision[x, 30], SetPrecision[
    y, 30]]}, {x, xRangeRC}, {y, yRangeRC}], 1];

rcNormal = Select[rcNormal, Element[#[[3]], Reals]&];

Export["../data/f32/wolfram/elliprc_data.csv", N[rcNormal], "CSV"];

(* PV domain: x >= 0, y < 0 *)

yRangeRCPV = Subdivide[-2, -muRC, 20];

rcPV = Flatten[Table[{x, y, CarlsonRC[SetPrecision[x, 30], SetPrecision[
    y, 30]]}, {x, xRangeRC}, {y, yRangeRCPV}], 1];

rcPV = Select[rcPV, Element[#[[3]], Reals]&];

Export["../data/f32/wolfram/elliprc_pv.csv", N[rcPV], "CSV"];

(* Carlson RD - Symmetric elliptic integral of the second kind *)

(* Normal domain: x >= 0, y >= 0, z > 0, at most one of x,y can be zero *)

xRangeRD = Union[Subdivide[muRD, 2, 9], {0}];

yRangeRD = Union[Subdivide[muRD, 2, 9], {0}];

zRangeRD = Subdivide[muRD, 2, 10];

validRD =
    Select[
        Flatten[Table[{x, y, z}, {x, xRangeRD}, {y, yRangeRD}, {z, zRangeRD
            }], 2]
        ,
        With[
            {x = #[[1]], y = #[[2]], z = #[[3]]}
            ,
            (* At most one of x, y can be zero *)
            !(x == 0 && y == 0)
        ]&
    ];

rdNormal =
    Select[
        Table[
            With[{x = data[[1]], y = data[[2]], z = data[[3]]},
                {x, y, z, CarlsonRD[SetPrecision[x, 30], SetPrecision[
                    y, 30], SetPrecision[z, 30]]}
            ]
            ,
            {data, validRD}
        ]
        ,
        Element[#[[4]], Reals]&
    ];

Export["../data/f32/wolfram/elliprd_data.csv", N[rdNormal], "CSV"];

(* Carlson RG - Symmetric elliptic integral of the second kind *)

(* Normal domain: x >= 0, y >= 0, z >= 0, at most one can be zero *)

xRangeRG = Union[Subdivide[muRG, 2, 9], {0}];

yRangeRG = Union[Subdivide[muRG, 2, 9], {0}];

zRangeRG = Union[Subdivide[muRG, 2, 9], {0}];

validRG =
    Select[
        Flatten[Table[{x, y, z}, {x, xRangeRG}, {y, yRangeRG}, {z, zRangeRG
            }], 2]
        ,
        With[
            {x = #[[1]], y = #[[2]], z = #[[3]]}
            ,
            (* At most one can be zero *)
            !(x == 0 && y == 0) && !(x == 0 && z == 0) && !(y == 0 &&
                 z == 0)
        ]&
    ];

rgNormal =
    Select[
        Table[
            With[{x = data[[1]], y = data[[2]], z = data[[3]]},
                {x, y, z, CarlsonRG[SetPrecision[x, 30], SetPrecision[
                    y, 30], SetPrecision[z, 30]]}
            ]
            ,
            {data, validRG}
        ]
        ,
        Element[#[[4]], Reals]&
    ];

Export["../data/f32/wolfram/elliprg_data.csv", N[rgNormal], "CSV"];

(* Carlson RF - Symmetric elliptic integral of the first kind *)

(* Normal domain: x >= 0, y >= 0, z >= 0, at most one can be zero *)

xRangeRF = Union[Subdivide[muRF, 2, 9], {0}];

yRangeRF = Union[Subdivide[muRF, 2, 9], {0}];

zRangeRF = Union[Subdivide[muRF, 2, 9], {0}];

validRF =
    Select[
        Flatten[Table[{x, y, z}, {x, xRangeRF}, {y, yRangeRF}, {z, zRangeRF
            }], 2]
        ,
        With[
            {x = #[[1]], y = #[[2]], z = #[[3]]}
            ,
            (* At most one can be zero *)
            !(x == 0 && y == 0) && !(x == 0 && z == 0) && !(y == 0 &&
                 z == 0)
        ]&
    ];

rfNormal =
    Select[
        Table[
            With[{x = data[[1]], y = data[[2]], z = data[[3]]},
                {x, y, z, CarlsonRF[SetPrecision[x, 30], SetPrecision[
                    y, 30], SetPrecision[z, 30]]}
            ]
            ,
            {data, validRF}
        ]
        ,
        Element[#[[4]], Reals]&
    ];

Export["../data/f32/wolfram/elliprf_data.csv", N[rfNormal], "CSV"];

(* Carlson RJ - Symmetric elliptic integral of the third kind *)

(* Normal domain: x >= 0, y >= 0, z >= 0, p > 0, at most one of x,y,z can be zero *)

xRangeRJ = Union[Subdivide[muRJ, 2, 5], {0}];

yRangeRJ = Union[Subdivide[muRJ, 2, 5], {0}];

zRangeRJ = Union[Subdivide[muRJ, 2, 5], {0}];

pRangeRJ = Subdivide[muRJ, 2, 5];

validRJ =
    Select[
        Flatten[Table[{x, y, z, p}, {x, xRangeRJ}, {y, yRangeRJ}, {z,
             zRangeRJ}, {p, pRangeRJ}], 3]
        ,
        With[
            {x = #[[1]], y = #[[2]], z = #[[3]], p = #[[4]]}
            ,
            (* At most one of x, y, z can be zero *)
            !(x == 0 && y == 0) && !(x == 0 && z == 0) && !(y == 0 &&
                 z == 0)
        ]&
    ];

rjNormal =
    Select[
        Table[
            With[{x = data[[1]], y = data[[2]], z = data[[3]], p = data
                [[4]]},
                {x, y, z, p, CarlsonRJ[SetPrecision[x, 30], SetPrecision[
                    y, 30], SetPrecision[z, 30], SetPrecision[p, 30]]}
            ]
            ,
            {data, validRJ}
        ]
        ,
        Element[#[[5]], Reals]&
    ];

Export["../data/f32/wolfram/elliprj_data.csv", N[rjNormal], "CSV"];

(* PV domain: x >= 0, y >= 0, z >= 0, p < 0, at most one of x,y,z can be zero *)

xRangeRJPV = Union[Subdivide[muRJPV, 2, 5], {0}];

yRangeRJPV = Union[Subdivide[muRJPV, 2, 5], {0}];

zRangeRJPV = Union[Subdivide[muRJPV, 2, 5], {0}];

pRangeRJPV = Subdivide[-2, -muRJPV, 5];

validRJPV =
    Select[
        Flatten[Table[{x, y, z, p}, {x, xRangeRJPV}, {y, yRangeRJPV},
             {z, zRangeRJPV}, {p, pRangeRJPV}], 3]
        ,
        With[
            {x = #[[1]], y = #[[2]], z = #[[3]], p = #[[4]]}
            ,
            (* At most one of x, y, z can be zero *)
            !(x == 0 && y == 0) && !(x == 0 && z == 0) && !(y == 0 &&
                 z == 0)
        ]&
    ];

rjPV =
    Select[
        Table[
            With[{x = data[[1]], y = data[[2]], z = data[[3]], p = data
                [[4]]},
                {x, y, z, p, CarlsonRJ[SetPrecision[x, 40], SetPrecision[
                    y, 40], SetPrecision[z, 40], SetPrecision[p, 40]]}
            ]
            ,
            {data, validRJPV}
        ]
        ,
        Element[#[[5]], Reals]&
    ];

Export["../data/f32/wolfram/elliprj_pv.csv", N[rjPV], "CSV"];

xyzRJPV = Subdivide[muRJPVSmall, muRJPV, 3];

pRangeRJPV = Subdivide[-2, -muRJPVSmall, 5];

validRJPVSmall =
    Select[
        Flatten[Table[{x, y, z, p}, {x, xyzRJPV}, {y, xyzRJPV}, {z, xyzRJPV
            }, {p, pRangeRJPV}], 3]
        ,
        With[
            {x = #[[1]], y = #[[2]], z = #[[3]], p = #[[4]]}
            ,
            (* At most one of x, y, z can be zero *)
            !(x == 0 && y == 0) && !(x == 0 && z == 0) && !(y == 0 &&
                 z == 0)
        ]&
    ];

rjPVSmall =
    Select[
        Table[
            With[{x = data[[1]], y = data[[2]], z = data[[3]], p = data
                [[4]]},
                {x, y, z, p, CarlsonRJ[SetPrecision[x, 40], SetPrecision[
                    y, 40], SetPrecision[z, 40], SetPrecision[p, 40]]}
            ]
            ,
            {data, validRJPVSmall}
        ]
        ,
        Element[#[[5]], Reals]&
    ];

Export["../data/f32/wolfram/elliprj_pv_small.csv", N[rjPVSmall], "CSV"
    ];
