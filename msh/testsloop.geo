
c1 = 0.2;
c2 = 0.05;
c3 = 0.01;

Point(1) = {0.0, 0.0, 0.0, c1};
Point(2) = {6.0, 0.0, 0, c1};
Point(3) = {6.0, 2.5, 0, c1};
Point(4) = {4, 2.5, 0, c3};
Point(5) = {3.5, 2.5, 0, c3};
Point(6) = {2.5, 2.5, 0, c3};
Point(7) = {2.0, 2.5, 0, c3};
Point(8) = {0.0, 2.5, 0, c1};
Point(9) = {5.9, 0.01, 0, c1};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {9, 2};




Curve Loop(1) = {1,2,3,4,5,6,7,8};


Plane Surface(1) = {1};
Curve{9} In Surface{1};


Physical Curve("Γᵍ₁") = {2,1,8};
Physical Curve("Γᵍ₂") = {5};
Physical Surface("Ω") = {1};
Physical Curve("Γᶜ") = {9};


Mesh.Algorithm = 1;
Mesh.MshFileVersion = 3;
Mesh 2;