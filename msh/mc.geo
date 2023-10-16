
c1 = 1;
c2 = 0.2;
c3 = 0.2;

Point(1) = {0.0, 0.0, 0.0, c2};
Point(2) = {13.3, 0.0, 0, c3};
Point(3) = {15.0, 0.0, 0, c1};
Point(4) = {15.0, 1.7, 0, c3};
Point(5) = {15, 50.0, 0, c1};
Point(6) = {0.0, 50.0, 0, c1};
Point(7) = {15, 20, 0, c1};
Point(8) = {0.0, 20.0, 0, c2};


Line(1) = {1, 2};
Circle(2) = {2, 3, 4};
Line(3) = {4, 7};
Line(4) = {7, 5};
Line(5) = {5, 6};
Line(6) = {6, 8};
Line(7) = {8, 1};




Curve Loop(1) = {1,2,3,4,5,6,7};


Plane Surface(1) = {1};


Curve{6} In Surface{1};

Physical Curve("Γᵗ") = {7,6};
Physical Curve("Γ²") = {1};
Physical Curve("Γ¹") = {3,4};
Physical Surface("Ω") = {1};


Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;