c1 = 0.003;
c2 = 0.0006;

Point(1) = {0, 0, 0, c1};
Point(2) = {0.08, 0, 0, c1};
Point(3) = {0.08, 0.16, 0, c1};
Point(4) = {0, 0.16, 0, c1};


//Point(7) = {0.035, 0.07625, 0,c2};
//Point(8) = {0.0465, 0.08375, 0,c2};
Point(7) = {0.036, 0.073072, 0,c2};
Point(8) = {0.044, 0.086928, 0,c2};
Point(9) = {0.0245, 0.065566, 0,c2};
Point(10) = {0.0555, 0.093560, 0,c2};
Point(11) = {0.001, 0.05, 0,c2};
Point(12) = {0.079, 0.11, 0,c2};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Line(5) = {7,9};
Line(6) = {8,10};
Line(7) = {8,7};
Line(8) = {9,11};
Line(9) = {10,12};


Curve Loop(1) = {1,2,3,4};


Plane Surface(1) = {1};
Curve{5,6,7,8,9} In Surface{1};

Physical Curve("Γᵍ₁") = {1};
Physical Curve("Γᵍ₂") = {3};
Physical Point("Γᵍ₃") = {1};
Physical Curve("Γᶜ") = {5,6};
Physical Surface("Ω") = {1};


// 生成网格
Mesh.Algorithm = 8;
Mesh.MshFileVersion = 3;
Mesh 2;
