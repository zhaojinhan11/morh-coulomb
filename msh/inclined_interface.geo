
c1 = 0.1;
c2 = 0.01;

Point(1) = {0, 0, 0, c1};
Point(2) = {1, 0, 0, c1};
Point(3) = {1, 1, 0, c1};
Point(4) = {0, 1, 0, c1};
Point(5) = {0.3, 0.33, 0, c2};
Point(6) = {0.7, 0.68, 0, c2};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Curve{5} In Surface{1};

Physical Curve("Γᵍ₁") = {1};
Physical Curve("Γᵍ₂") = {3};
Physical Curve("Γᶜ") = {5};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;