// Gmsh project created on Mon Jul 10 23:03:10 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0.0, 0.0, 0.0, 1};
//+
Point(2) = {13.6, 0.0, 0, 0.1};
//+
Point(3) = {15.0, 0.0, 0, 1};
//+
Point(4) = {15.0, 1.4, 0, 1};
//+
Point(5) = {15, 50.0, 0, 1};
//+
Point(6) = {0.0, 50.0, 0, 1};
//+
Line(1) = {1, 2};

//+
Circle(2) = {2, 3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 1};
//+
Transfinite Curve {3,  5} = 50 Using Progression 1;
//+
Transfinite Curve {4, 1，2} = 100 Using Progression 1;
//+
Transfinite Curve {2} = 50 Using Progression 1;
//+
Curve Loop(1) = {1，2，3，4，5};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1};
//+
Physical Curve("Γᵗ") = {5};
//+
Physical Curve("Γ²") = {1};
//+
Physical Curve("Γ¹") = {3};
//+
Physical Surface("Ω") = {1};
