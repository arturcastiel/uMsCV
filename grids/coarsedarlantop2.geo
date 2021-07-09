cl__1 = 3;
h = 0.49;
//h = 0.5;
H = h/2;
// 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
Mesh.Algorithm = 2;
Mesh.CharacteristicLengthFactor = 2;
diva = 5;
divb = 4;
Point(1) = {0, 0, 0, cl__1};
Point(2) = {1, 0, 0, cl__1};
Point(3) = {1, 1, 0, cl__1};
Point(4) = {0, 1, 0, cl__1};
Point(5) = {0.5 - H, 0.5 - H, 0, cl__1};
Point(6) = {0.5+H, 0.5-H, 0, cl__1};
Point(7) = {0.5+H, 0.5+H, 0, cl__1};
Point(8) = {0.5-H, 0.5+H, 0, cl__1};
Line(1) = {1, 2};

Line(2) = {2, 3};

Line(3) = {3, 4};

Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};

Plane Surface(2) = {2};

Transfinite Line {1} = diva Using Progression 1.05;
Transfinite Line {2} = diva Using Progression 0.98;
Transfinite Line {3} = diva Using Progression 1;
Transfinite Line {4} = diva Using Progression 1;


Transfinite Line {5} = divb Using Progression 1;
Transfinite Line {6} = divb Using Progression 1;
Transfinite Line {7} = divb Using Progression 1;
Transfinite Line {8} = divb Using Progression 1;

//Line {5} In Surface {2};
//Line {6} In Surface {2};
//Line {7} In Surface {2};
//Line {8} In Surface {2};
Recombine Surface {1};
Recombine Surface {2};

Physical Point(201) = {1, 2, 3, 4};
Physical Line(201) = {1, 2, 3, 4};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
