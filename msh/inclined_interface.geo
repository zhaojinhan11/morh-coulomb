c1 = 0.013;
c2 = 0.013;

Point(1) = {0, 0, 0, c1};
Point(2) = {1, 0, 0, c1};
Point(3) = {1, 1, 0, c1};
Point(4) = {0, 1, 0, c1};
Point(11) = {0.05, 0, 0, c1};
Point(12) = {0.1, 0, 0, c1};
Point(13) = {0.15, 0, 0, c1};
Point(14) = {0.2, 0, 0, c1};
Point(31) = {1, 0.95, 0, c1};


Point(7) = {0.3, 0.3, 0,c2};
Point(8) = {0.7, 0.7, 0,c2};
Point(9) = {0.000, 0.000, 0,c2};
Point(10) = {0.999,0.999, 0,c2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Line(5) = {7,8};
Line(6) = {7,9};
Line(7) = {8,10};



Curve Loop(1) = {1,2,3,4};


Plane Surface(1) = {1};


Curve{7} In Surface{1};

Physical Curve("Γᵍ₁") = {1};
Physical Curve("Γᵍ₂") = {3};
Physical Curve("Γᶜ") = {5};
Physical Surface("Ω") = {1};


// 生成网格
Mesh.Algorithm = 1;
Mesh.MshFileVersion = 3;
Mesh 2;
