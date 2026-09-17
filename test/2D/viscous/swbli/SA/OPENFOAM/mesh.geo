// =====================================================
// SWBLI – Shock/Boundary Layer Interaction (OpenFOAM)
// Extruded (1-layer) version of ../script.geo for
// gmshToFoam. Differences vs the MOSE script:
//   * top wall split at x = 0 (symTop upstream / topWall
//     downstream) to replicate the MOSE [up] patch split
//   * extrusion of one hex layer in z + Physical groups
// Keep the node counts below in sync with ../script.geo.
// =====================================================

f = 1.41;

bl_h_1 = 0.0025;
bl_h_2 = 0.008;

h_top = 0.004/f;

rate_1 = 1.1 ^ (1.0/f);
rate_2 = 1.05 ^ (1.0/f);

n_bl_1 = 73;   // 72 cells = 9×8
n_bl_2 = 97;   // 96 cells = 12×8
n_in   = 89;   // 88 cells = 11×8
n_sh   = 225;  // 224 cells = 28×8
n_out  = 89;   // 88 cells = 11×8

dz = 0.01;   // extrusion thickness (arbitrary: front/back are empty)

// -------------------- POINTS --------------------
// Row 1 – bottom wall (y = 0)
Point(1)  = {-0.01,   0, 0, h_top/1.5};
Point(2)  = {0.0,     0, 0, h_top/1.5};
Point(3)  = {0.16544, 0, 0, h_top/1.5};
Point(4)  = {0.31844, 0, 0, h_top/5};
Point(5)  = {0.37,    0, 0, h_top/5};
Point(6)  = {0.523,   0, 0, h_top};

// Row 2 – bottom BL top (non-uniform height)
Point(11) = {-0.01,   3.5*bl_h_1,  0, h_top/1.5};
Point(12) = {0.0,     3.5*bl_h_1,  0, h_top/1.5};
Point(13) = {0.16544, 2.25*bl_h_1, 0, h_top/1.5};
Point(14) = {0.31844, bl_h_1,      0, h_top/5};
Point(15) = {0.37,    bl_h_1,      0, h_top/5};
Point(16) = {0.523,   bl_h_1,      0, h_top};

// Row 3 – top BL bottom
Point(21) = {-0.01,   0.115    - bl_h_2, 0, h_top};
Point(22) = {0.023,   0.115    - bl_h_2, 0, h_top/4};
Point(23) = {0.31844, 0.062906 - bl_h_2, 0, h_top};
Point(24) = {0.523,   0.062906 - bl_h_2, 0, h_top};
Point(25) = {0.37,    0.062906 - bl_h_2, 0, h_top};

// Row 4 – top wall (split at x = 0: symmetry / isothermal)
Point(31) = {-0.01,   0.115,    0, h_top};
Point(36) = {0.0,     0.115,    0, h_top/4};
Point(32) = {0.023,   0.115,    0, h_top/4};
Point(33) = {0.31844, 0.062906, 0, h_top};
Point(34) = {0.523,   0.062906, 0, h_top};
Point(35) = {0.37,    0.062906, 0, h_top};

// -------------------- LINES --------------------
// Bottom wall (left → right)
Line(1)  = {1,  2};
Line(2)  = {2,  3};
Line(3)  = {3,  4};
Line(4)  = {4,  5};
Line(5)  = {5,  6};

// Bottom BL top (left → right)
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};

// Top BL bottom (left → right)
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 25};
Line(24) = {25, 24};

// Top wall (left → right); 31 = symmetry part, 131 = wall part
Line(31)  = {31, 36};
Line(131) = {36, 32};
Line(32)  = {32, 33};
Line(33)  = {33, 35};
Line(34)  = {35, 34};

// Vertical lines
Line(41) = {1,  11};   // inlet,  bottom BL
Line(42) = {6,  16};   // outlet, bottom BL
Line(43) = {11, 21};   // inlet,  interior
Line(44) = {16, 24};   // outlet, interior
Line(45) = {21, 31};   // inlet,  top BL
Line(46) = {24, 34};   // outlet, top BL
Line(51) = {14, 23};   // interior split, x=0.31844
Line(52) = {15, 25};   // interior split, x=0.37000

// -------------------- SURFACES --------------------
// Bottom boundary layer
Curve Loop(1) = {1, 2, 3, 4, 5, 42, -15, -14, -13, -12, -11, -41};
Plane Surface(1) = {1};

// Interior – Block A (x = -0.01 → 0.31844)
Curve Loop(3) = {11, 12, 13, 51, -22, -21, -43};
Plane Surface(3) = {3};

// Interior – Block B (x = 0.31844 → 0.370)
Curve Loop(4) = {51, 23, -52, -14};
Plane Surface(4) = {4};

// Interior – Block C (x = 0.37000 → 0.523)
Curve Loop(5) = {52, 24, -44, -15};
Plane Surface(5) = {5};

// Top boundary layer
Curve Loop(2) = {21, 22, 23, 24, 46, -34, -33, -32, -131, -31, -45};
Plane Surface(2) = {2};

// -------------------- TRANSFINITE --------------------
// --- Horizontal curves shared by BL and interior rows ---
Transfinite Curve {1,  11}         = n_in;              // x: -0.01 → 0.0
Transfinite Curve {2,  3,  12, 13} = n_sh;              // x: 0.0 → 0.16544 → 0.31844
Transfinite Curve {4,  14,  23, 33} = n_sh Using Bump 6;// x: 0.31844 → 0.37 (shock)
Transfinite Curve {5,  15,  24, 34} = n_out;            // x: 0.37 → 0.523

// Top BL + interior top boundary
Transfinite Curve {21} = n_in;
Transfinite Curve {22, 32} = 449;   // = 2*n_sh - 1, matches Block-A bottom

// Top wall split at x=0: the two pieces must sum back to n_in
// (n_l + n_r - 1 = n_in) so the composite side still matches curve 21.
// Split proportionally to length: x -0.01→0 is 10/33 of the -0.01→0.023 span.
n_in_l = Round( (n_in - 1) * 10.0/33.0 ) + 1;   // x: -0.01 → 0.0
n_in_r = n_in - n_in_l + 1;                     // x:  0.0  → 0.023
Transfinite Curve {31}  = n_in_l;
Transfinite Curve {131} = n_in_r;

// === BOTTOM BL (Surface 1) ===
Transfinite Curve {41, 42}  = n_bl_1 Using Progression rate_1;
Transfinite Surface {1} = {1, 6, 16, 11};

// === TOP BL (Surface 2) ===
Transfinite Curve {-45, -46} = n_bl_2 Using Progression rate_2;
Transfinite Surface {2} = {21, 24, 34, 31};

// === INTERIOR Blocks A/B/C ===
Transfinite Curve {43, 51, 52, 44} = 113 Using Bump 0.25;
Transfinite Surface {3} = {11, 14, 23, 21};
Transfinite Surface {4} = {14, 15, 25, 23};
Transfinite Surface {5} = {15, 16, 24, 25};

Recombine Surface {1, 2, 3, 4, 5};

// -------------------- EXTRUSION --------------------
// out[0] = back-plane copy of the surface, out[1] = volume,
// out[2...] = lateral surfaces, one per curve, in loop order.
o1[] = Extrude {0, 0, dz} { Surface{1}; Layers{1}; Recombine; };
o3[] = Extrude {0, 0, dz} { Surface{3}; Layers{1}; Recombine; };
o4[] = Extrude {0, 0, dz} { Surface{4}; Layers{1}; Recombine; };
o5[] = Extrude {0, 0, dz} { Surface{5}; Layers{1}; Recombine; };
o2[] = Extrude {0, 0, dz} { Surface{2}; Layers{1}; Recombine; };

// -------------------- PHYSICAL GROUPS --------------------
// Surface 1 loop: {1,2,3,4,5,42,-15,-14,-13,-12,-11,-41}
//   o1[2..6] = bottom wall, o1[7] = outlet, o1[13] = inlet
// Surface 3 loop: {11,12,13,51,-22,-21,-43}
//   o3[8] = inlet
// Surface 5 loop: {52,24,-44,-15}
//   o5[4] = outlet
// Surface 2 loop: {21,22,23,24,46,-34,-33,-32,-131,-31,-45}
//   o2[6] = outlet, o2[7..10] = top wall, o2[11] = symTop, o2[12] = inlet

Physical Surface("inlet")        = {o1[13], o3[8], o2[12]};
Physical Surface("outlet")       = {o1[7],  o5[4], o2[6]};
Physical Surface("bottomWall")   = {o1[2], o1[3], o1[4], o1[5], o1[6]};
Physical Surface("topWall")      = {o2[7], o2[8], o2[9], o2[10]};
Physical Surface("symTop")       = {o2[11]};
Physical Surface("frontAndBack") = {1, 2, 3, 4, 5, o1[0], o2[0], o3[0], o4[0], o5[0]};

Physical Volume("fluid") = {o1[1], o2[1], o3[1], o4[1], o5[1]};

// gmshToFoam requires the legacy MSH2 format
Mesh.MshFileVersion = 2.2;
