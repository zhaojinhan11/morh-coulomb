
c1 = 0.03;
c2 = 0.01;

Point(1) = {0, 0, 0, c1};
Point(2) = {1, 0, 0, c1};
Point(3) = {1, 1, 0, c1};
Point(4) = {0, 1, 0, c1};
Point(5) = {0, 0.5, 0, c1};
Point(6) = {0.5, 0.5, 0, c2};
Point(7) = {1, 0.5, 0, c1};
Point(8) = {0.499, 0.499, 0, c2};
Point(9) = {0.9, 0.01, 0, c2};




Line(1) = {1,2};
Line(2) = {2,7};
Line(3) = {7,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {5,1};
Line(9) = {8,9};



Curve Loop(1) = {1,2,-7,-6,8};
Curve Loop(2) = {3,4,5,6,7};


Plane Surface(1) = {1};
Plane Surface(2) = {2};
Curve{9} In Surface{1};

Physical Curve("Γᵍ₁") = {1};
Physical Curve("Γᵍ₂") = {4};
Physical Curve("Γᶜ") = {6};
Physical Surface("Ω") = {1,2};

Transfinite Curve{6}=70;
Transfinite Curve{9}=300;
Transfinite Surface{1};


Mesh.Algorithm = 8;
Mesh.MshFileVersion = 3;
Mesh 2;
