
c1 = 1;
c2 = 0.1;

Point(1) = {0, 0, 0, c2};
Point(2) = {20, 0, 0, c1};
Point(3) = {20, 10, 0, c1};
Point(4) = {14, 10, 0, c2};
Point(5) = {10, 10, 0, c1};
//Point(6) = {19.9, 0.1, 0, c2};
//Point(7) = {19.9, 0.2, 0, c2};



Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};
Line(6) = {4,1};
//Line(7) = {6,7};



Curve Loop(1) = {1,2,3,4,5};



Plane Surface(1) = {1};


Curve{6} In Surface{1};

Physical Curve("Γᵍ₁") = {1};
Physical Curve("Γᵍ₂") = {4};
Physical Curve("Γᵍ₃") = {2};
//Physical Curve("Γᶜ") = {7};
Physical Surface("Ω") = {1};


Mesh.Algorithm = 1;
Mesh.MshFileVersion = 3;
Mesh 2;