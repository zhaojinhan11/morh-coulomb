
a = 48.0;
b = 12.0;
n = 4;

Point(1) = {0.0, -b/2, 0.0};
Point(2) = {  a, -b/2, 0.0};
Point(3) = {  a,  b/2, 0.0};
Point(4) = {0.0,  b/2, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Transfinite Curve{1,3} = 4*n+1;
Transfinite Curve{2,4} = n+1;

Physical Curve("Γᵗ") = {2};
Physical Curve("Γᵍ") = {4};
Physical Surface("Ω") = {1};

Transfinite Surface{1};

Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
RecombineMesh;


