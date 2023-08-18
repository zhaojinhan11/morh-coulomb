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
Line(6) = {1, 5};
//+
Line(7) = {6, 3};
//+
Transfinite Curve {5, 4, 3, 2, 6, 7} = 30 Using Progression 1;
//+
Transfinite Curve {1} = 130 Using Progression 1;
//+
Curve Loop(1) = {4, 5, 6, 1, 7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -7, -1, -6, 2};
//+
Plane Surface(2) = {2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Line{1} In Surface{1};
//+
Line{1} In Surface{2};
//+
Physical Curve("Γ") = {2};
//+
Physical Curve("Γᵛ") = {1};
//+
Physical Curve("Γᵍ") = {4};
//+
Physical Surface("Ω") = {1,2};
