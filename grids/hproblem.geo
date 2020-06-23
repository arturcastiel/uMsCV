cl__1 = 1e+22;
ep = 43;
h = 0.5;
x = 8;
x = x + 1;
// 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
Mesh.Algorithm =4;
Mesh.CharacteristicLengthFactor = 8;
qi = 0.5 - (h/2);
qf = 0.5 + (h/2);
Point(1) = {0, 0, 0, ep};
Point(2) = {1, 0, 0, ep};
Point(3) = {1, 1, 0, ep};
Point(4) = {0, 1, 0, ep};
Point(5) = {qi, qi, 0, ep};
Point(6) = {qf, qi, 0, ep};
Point(7) = {qf, qf, 0, ep};
Point(8) = {qi, qf, 0, ep};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line Loop(1) = {1, 2, 3, 4, -8, -7, -6, -5};
Plane Surface(1) = {1};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Line {5} In Surface {2};
Line {6} In Surface {2};
Line {7} In Surface {2};
Line {8} In Surface {2};




//Transfinite Line{1} = 7Using Progression 1.15;
//Transfinite Line{4} = 7Using Progression 0.8696;


//Transfinite Line{1} = 6Using Progression 1.15;
//Transfinite Line{4} = 6Using Progression 0.8696;



//Transfinite Line{1} = 7Using Progression 1.15;
//Transfinite Line{4} = 7Using Progression 0.8696;

//Transfinite Line{2} = 6Using Progression 0.8696;
//Transfinite Line{3} = 6Using Progression 1.15;

//Transfinite Line{2} = 6Using Progression 0.8696;
//Transfinite Line{3} = 6Using Progression 1.15;

//Transfinite Line{1,2,3,4} = x Using Progression 1;

//Transfinite Line{2,3} = 6Using Progression 1.1;
//Transfinite Line{5,6,7,8} = (x ) Using Progression 1;
//Transfinite Surface{1} = {1};
//Transfinite Surface{1} = {1,2,3,4};
//Transfinite Surface{2} = {5,6,7,8};
//Transfinite Line{1,3} = 4Using Progression 1;
//Transfinite Line{2,4} = 4Using Progression 1;

// Recombine the triangles into quads:
//Transfinite Surface{1} = {1,2,3,4};
Recombine Surface{1};
Recombine Surface{2};
Mesh.Smoothing = 100;
Physical Point(101) = {1, 2, 3, 4};
Physical Line(101) = {1, 2, 3, 4};
Physical Surface(1) = {1};
Physical Surface(2) = {2};