cl__1 = 1e+22;
ep = 0;
h = 0.5;
x = 8;
x = x + 1;
// 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
Mesh.Algorithm =2;
//Mesh.CharacteristicLengthFactor = 5;
qi = 0.5 - (h/2);
qf = 0.5 + (h/2);
Point(1) = {0, 0, 0, ep};
Point(2) = {1, 0, 0, ep};
Point(3) = {1, 1, 0, ep};
Point(4) = {0, 1, 0, ep};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};





Transfinite Line{1} = 28 Using Progression 1;

Transfinite Line{2} =28 Using Progression 1;

Transfinite Line{3} = 29 Using Progression 1;
Transfinite Line{4} = 38 Using Progression 1;

// Recombine the triangles into quads:
//Transfinite Surface{1} = {1,2,3,4};
//Recombine Surface{1};
//Recombine Surface{2};
Mesh.Smoothing = 100;
Physical Point(201) = {1, 2, 3, 4};
Physical Line(201) = {1, 2, 3, 4};
Physical Surface(1) = {1};
//Physical Surface(2) = {2};
