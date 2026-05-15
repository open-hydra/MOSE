// =====================================================
// SWBLI – Shock/Boundary Layer Interaction (Geometry)
// =====================================================

f = 1.41;

bl_h_1 = 0.0025;
bl_h_2 = 0.008;

h_top = 0.004/f;

n_bl_1 = 50*f+1;
n_bl_2 = 31*f+1;
n_in   = 60*f+1;
n_sh   = 160*f+1;
n_out  = 64*f+1;
rate_1 = 1.1 ^ (1.0/f);
rate_2 = 1.25 ^ (1.0/f);

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

// Row 4 – top wall
Point(31) = {-0.01,   0.115,    0, h_top};
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

// Top wall (left → right)
Line(31) = {31, 32};
Line(32) = {32, 33};
Line(33) = {33, 35};
Line(34) = {35, 34};

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
Curve Loop(2) = {21, 22, 23, 24, 46, -34, -33, -32, -31, -45};
Plane Surface(2) = {2};

// -------------------- TRANSFINITE --------------------
// --- Horizontal curves shared by BL and interior rows ---
// (same count per x-segment ensures all opposite sides match)
Transfinite Curve {1,  11}        = n_in;             // x: -0.01 → 0.0
Transfinite Curve {2,  3,  12, 13} = n_sh;            // x: 0.0 → 0.16544 → 0.31844
Transfinite Curve {4,  14,  23, 33}    = n_sh Using Bump 6;// x: 0.31844 → 0.37 (shock)
Transfinite Curve {5,  15,  24, 34}   = n_out;            // x: 0.37 → 0.523

// Top BL + interior top boundary
// N21=n_in, N22=2*n_sh-1, N23=n_sh+n_out-1
// ensures Block-A top = n_in+2*n_sh-2 = Block-A bottom ✓
//         Block-B top = n_sh+n_out-1  = Block-B bottom ✓
//         Top-BL bottom total = top-wall total           ✓
Transfinite Curve {21, 31} = n_in;
Transfinite Curve {22, 32} = 2*n_sh - 2;

// === BOTTOM BL (Surface 1) ===
Transfinite Curve {41, 42}  = n_bl_1 Using Progression rate_1;
Transfinite Surface {1} = {1, 6, 16, 11};

// === TOP BL (Surface 2) ===
Transfinite Curve {-45, -46} = n_bl_2 Using Progression rate_2;
Transfinite Surface {2} = {21, 24, 34, 31};

// === INTERIOR Block A (Surface 3, corners 11-14-23-21) ===
// === INTERIOR Block B (Surface 4, corners 14-16-24-23) ===
Transfinite Curve {43, 51, 52, 44} = 1.5*(0.115 - bl_h_1 - bl_h_2) / h_top Using Bump 0.25;
Transfinite Surface {3} = {11, 14, 23, 21};
Transfinite Surface {4} = {14, 15, 25, 23};
Transfinite Surface {5} = {15, 16, 24, 25};

Recombine Surface {1, 2, 3, 4, 5};

