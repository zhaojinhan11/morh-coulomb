
c1 = 0.013;
c2 = 0.001;

Point(1) = {0, 0, 0, c1};
Point(2) = {0.8, 0, 0, c1};
Point(3) = {0.8, 0.2, 0, c1};
Point(4) = {0, 0.2, 0, c1};
Point(5) = {0.4, 0.04, 0, c2};
Point(6) = {0.4, 0.16, 0, c2};
Point(7) = {0.4, 0, 0, c2};
Point(8) = {0.4, 0.2, 0, c2};
Point(9) = {0.4, 0.15, 0, c2};
Point(10) = {0.32, 0.05, 0, c2};



Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {8,4};
Line(6) = {4,1};
Line(7) = {7,5};
Line(8) = {6,8};

Curve Loop(1) = {1,2,3,4};


Plane Surface(1) = {1};
Curve{7,8} In Surface{1};

Physical Curve("Γᵗ") = {2};
Physical Curve("Γᵍ") = {4};
Physical Surface("Ω") = {1};

Transfinite Surface{1};


Mesh.Algorithm = 4;
Mesh.MshFileVersion = 5;
Mesh 4;


