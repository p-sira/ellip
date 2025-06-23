(*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 *)

epsilon = 2.22*^-16;

(* Mu helps prevent answers approaching infinity, making the comparison fairer. *)

muRC = epsilon;

muRD = 50 * epsilon;

muRG = 50 * epsilon;

muRF = epsilon;

muRJ = 100 * epsilon;

(* Carlson RC - Degenerate symmetric elliptic integral of RF *)

(* Normal domain: x >= 0, y > 0 *)

xRangeRC = Union[Subdivide[muRC, 2, 30], {0}];

yRangeRC = Subdivide[muRC, 2, 30];

rcNormal = Flatten[Table[{x, y, CarlsonRC[x, y]}, {x, xRangeRC}, {y, 
  yRangeRC}], 1];

rcNormal = Select[rcNormal, Element[#[[3]], Reals]&];

Export["../data/wolfram/elliprc_data.csv", N[rcNormal], "CSV"];

(* PV domain: x >= 0, y < 0 *)

yRangeRCPV = Subdivide[-2, -muRC, 20];

rcPV = Flatten[Table[{x, y, CarlsonRC[x, y]}, {x, xRangeRC}, {y, yRangeRCPV
  }], 1];

rcPV = Select[rcPV, Element[#[[3]], Reals]&];

Export["../data/wolfram/elliprc_pv.csv", N[rcPV], "CSV"];

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
        {x, y, z, CarlsonRD[x, y, z]}
      ]
      ,
      {data, validRD}
    ]
    ,
    Element[#[[4]], Reals]&
  ];

Export["../data/wolfram/elliprd_data.csv", N[rdNormal], "CSV"];

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
      !(x == 0 && y == 0) && !(x == 0 && z == 0) && !(y == 0 && z == 
        0)
    ]&
  ];

rgNormal =
  Select[
    Table[
      With[{x = data[[1]], y = data[[2]], z = data[[3]]},
        {x, y, z, CarlsonRG[x, y, z]}
      ]
      ,
      {data, validRG}
    ]
    ,
    Element[#[[4]], Reals]&
  ];

Export["../data/wolfram/elliprg_data.csv", N[rgNormal], "CSV"];

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
      !(x == 0 && y == 0) && !(x == 0 && z == 0) && !(y == 0 && z == 
        0)
    ]&
  ];

rfNormal =
  Select[
    Table[
      With[{x = data[[1]], y = data[[2]], z = data[[3]]},
        {x, y, z, CarlsonRF[x, y, z]}
      ]
      ,
      {data, validRF}
    ]
    ,
    Element[#[[4]], Reals]&
  ];

Export["../data/wolfram/elliprf_data.csv", N[rfNormal], "CSV"];

(* Carlson RJ - Symmetric elliptic integral of the third kind *)

(* Normal domain: x >= 0, y >= 0, z >= 0, p > 0, at most one of x,y,z can be zero *)

xRangeRJ = Union[Subdivide[muRJ, 2, 5], {0}];

yRangeRJ = Union[Subdivide[muRJ, 2, 5], {0}];

zRangeRJ = Union[Subdivide[muRJ, 2, 5], {0}];

pRangeRJ = Subdivide[muRJ, 2, 5];

validRJ =
  Select[
    Flatten[Table[{x, y, z, p}, {x, xRangeRJ}, {y, yRangeRJ}, {z, zRangeRJ
      }, {p, pRangeRJ}], 3]
    ,
    With[
      {x = #[[1]], y = #[[2]], z = #[[3]], p = #[[4]]}
      ,
      (* At most one of x, y, z can be zero *)
      !(x == 0 && y == 0) && !(x == 0 && z == 0) && !(y == 0 && z == 
        0)
    ]&
  ];

rjNormal =
  Select[
    Table[
      With[{x = data[[1]], y = data[[2]], z = data[[3]], p = data[[4]]
        },
        {x, y, z, p, CarlsonRJ[x, y, z, p]}
      ]
      ,
      {data, validRJ}
    ]
    ,
    Element[#[[5]], Reals]&
  ];

Export["../data/wolfram/elliprj_data.csv", N[rjNormal], "CSV"];

(* PV domain: x >= 0, y >= 0, z >= 0, p < 0, at most one of x,y,z can be zero *)

pRangeRJPV = Subdivide[-2, -muRJ, 5];

validRJPV =
  Select[
    Flatten[Table[{x, y, z, p}, {x, xRangeRJ}, {y, yRangeRJ}, {z, zRangeRJ
      }, {p, pRangeRJPV}], 3]
    ,
    With[
      {x = #[[1]], y = #[[2]], z = #[[3]], p = #[[4]]}
      ,
      (* At most one of x, y, z can be zero *)
      !(x == 0 && y == 0) && !(x == 0 && z == 0) && !(y == 0 && z == 
        0)
    ]&
  ];

rjPV =
  Select[
    Table[
      With[{x = data[[1]], y = data[[2]], z = data[[3]], p = data[[4]]
        },
        {x, y, z, p, CarlsonRJ[x, y, z, p]}
      ]
      ,
      {data, validRJPV}
    ]
    ,
    Element[#[[5]], Reals]&
  ];

Export["../data/wolfram/elliprj_pv.csv", N[rjPV], "CSV"];
