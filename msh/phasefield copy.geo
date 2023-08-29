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

//+
Curve Loop(1) = {2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+

Line {1} In Surface{1};


//+
Physical Curve("Γ") = {2};
//+
Physical Curve("Γᵛ") = {1};
//+
Physical Curve("Γᵍ") = {4};
//+
Physical Surface("Ω") = {1};
