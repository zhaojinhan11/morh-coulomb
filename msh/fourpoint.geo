z = 20;
c1 = 0.01;
c2 = 0.01;


Point(1) = {0, 0, 0, c1};
Point(2) = {0.8, 0, 0, c1};
Point(3) = {0.8, 0.2, 0, c1};
Point(4) = {0, 0.2, 0, c1};
Point(5) = {0.4, 0.04, 0, c2};
Point(6) = {0.4, 0.16, 0, c2};
Point(7) = {0.4, 0, 0, c2};
Point(8) = {0.4, 0.2, 0, c2};
Point(9) = {0.399, 0.159, 0};
Point(10) = {0.329, 0.059, 0};
Point(11) = {0.471, 0.141, 0};
Point(12) = {0.401, 0.0399, 0, c2};
Point(13) = {0.35, 0, 0, c2};


Line(1) = {1,7};
Line(2) = {7,2};
Line(3) = {2,3};
Line(4) = {3,8};
Line(5) = {8,4};
Line(6) = {4,1};
Line(7) = {7,5};
Line(8) = {5,6};
Line(9) = {6,8};
Line(10) = {9,10};
Line(11) = {11,12};


Curve Loop(1) = {1,7,8,9,5,6};
Curve Loop(2) = {2,3,4,-9,-8,-7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Transfinite Curve{9,7} = 40;
Transfinite Curve{10,11} = 100;
Transfinite Surface{1};
Transfinite Surface{2};
Curve{10} In Surface{1};
Curve{11} In Surface{2};

Physical Point("Γᵍ₃") = {2};
Physical Curve("Γᵗ") = {2};
Physical Curve("Γᵍ") = {4};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 8;
Mesh.MshFileVersion = 3;
Mesh 2;


