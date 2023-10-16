
c1 = 0.05;
c2 = 0.005;

Point(1) = {0, 0, 0, c1};
Point(2) = {1, 0, 0, c1};
Point(3) = {1, 0.8, 0, c2};
Point(4) = {1, 1, 0, c1};
Point(5) = {0, 1, 0, c1};
Point(6) = {0, 0.2, 0, c2};
Point(7) = {0.3, 0.33, 0, c2};
Point(8) = {0.7, 0.68, 0, c2};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {7,8};
Line(8) = {6,7};
Line(9) = {8,3};




Curve Loop(1) = {1,2,3,4,5,6};


Plane Surface(1) = {1};


Curve{7,8,9} In Surface{1};

Physical Curve("Γᵍ₁") = {1};
Physical Curve("Γᵍ₂") = {4};
Physical Curve("Γᶜ") = {7};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;