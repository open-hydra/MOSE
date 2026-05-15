//+
d = 1.0;
r = 0.5*d;
t = 2;

Point(1) = {0, 0+t, 0, 1.0};
Point(2) = {0, -d/2+t, 0, 1.0};
Point(3) = {0, d/2+t, 0, 1.0};

Circle(1) = {3, 1, 2};

Point(4) = {0, 2.5*r+t, 0, 1.0};
Point(5) = {0, -2.5*r+t, 0, 1.0};

Line(2) = {3, 4};
Line(3) = {2, 5};
Circle(4) = {4, 1, 5};

Transfinite Curve {2, 3} = 41 Using Progression 1.0;
Transfinite Curve {1, 4} = 161 Using Progression 1.0;

Curve Loop(1) = {4, -3, -1, 2};
Plane Surface(1) = {1};
Transfinite Surface {1} = {2, 3, 4, 5};
