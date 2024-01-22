
c1 = 0.03;
c2 = 0.009;

Point(1) = {0, 0, 0, c1};
Point(2) = {1, 0, 0, c1};
Point(3) = {1, 1, 0, c1};
Point(4) = {0, 1, 0, c1};
Point(5) = {0, 0.4, 0, c1};
Point(6) = {1, 0.6, 0, c1};



Line(1) = {1,2};
Line(2) = {2,6};
Line(3) = {6,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,1};
Line(7) = {5,6};






Curve Loop(1) = {5,4,3,7};
Curve Loop(2) = {1,2,-7,6};


Plane Surface(1) = {1};
Plane Surface(2) = {2};

//Curve{7} In Surface{1};

Physical Curve("Γᵍ₁") = {1};
Physical Curve("Γᵍ₂") = {4};
Physical Curve("Γᵍ₃") = {3};
Physical Curve("Γᶜ") = {7};
Physical Surface("Ω") = {1,2};


Mesh.Algorithm = 1;
Mesh.MshFileVersion = 3;
Mesh 2;