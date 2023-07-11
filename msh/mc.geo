// Gmsh project created on Mon Jul 10 23:03:10 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {30, 0, 0, 1.0};
//+
Point(3) = {30, 100, 0, 1.0};
//+
Point(4) = {0, 100, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {15.0, 50.0, 0.0, 1.4, 0, 2*Pi};
//+
Transfinite Curve {4, 2, 5} = 50 Using Progression 1;
//+
Transfinite Curve {3, 1} = 15 Using Progression 1;
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
//+
Transfinite Surface {1};
//+
Physical Curve("Γᵗ") = {2};
//+
Physical Curve("Γ") = {4};
//+
Physical Surface("Ω") = {1};