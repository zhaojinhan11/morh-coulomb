// Gmsh project created on Fri Aug 04 15:56:41 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.05};
//+
Point(2) = {1, 0, 0, 0.05};
//+
Point(3) = {1, 1, 0, 0.05};
//+
Point(4) = {0, 1, 0, 0.05};
//+
Point(5) = {0.3, 0.33, 0, 0.01};
//+
Point(6) = {0.7, 0.68, 0, 0.01};
//+
Line(1) = {5, 6};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 1};
//+
Line(6) = {4, 6};
//+
Line(7) = {4, 5};
//+
Line(8) = {2, 5};
//+
Line(9) = {2, 6};
//+
Transfinite Curve {2，3，4，5，6，7，8，9} = 10 Using Progression 1;
//+
Transfinite Curve {1} = 30 Using Progression 1;
//+
Curve Loop(1) = {2, 8, 7, 5};
//+
Curve Loop(2) = {3, 4, 6, 9};
//+
Curve Loop(3) = {1, 6, 7};
//+
Curve Loop(4) = {1, 8, 9};
//+
Plane Surface(1) = {1};
//+
Plane Surface(2) = {2};
//+
Plane Surface(3) = {3};
//+
Plane Surface(4) = {4};
//
Physical Curve("Γ") = {2};
//+
Physical Curve("Γᵛ") = {1};
//+
Physical Curve("Γᵍ") = {4};

Transfinite Surface {3};
//+
Transfinite Surface {2};
//+
Transfinite Surface {4};
//+
Transfinite Surface {1};
//+
Physical Surface("Ω") = {1};
//+
Physical Surface("Ω") = {2};
//+
Physical Surface("Ω") = {3};
//+
Physical Surface("Ω") = {4};
