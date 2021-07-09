// Gmsh project created on Mon Jun 22 20:44:21 2020
cl1 = 0.045;
cl2 = 0.009;
cl3 = 0.009;
cl4 = 0.009;
cl2 = 0.004;
cl3 = 0.004;
cl4 = 0.004;
// 0.015
r1 = 0.015;
r2 = 0.015;
r3 = 0.015;
x1 = 0.3;
y1 = 0.3;
x2 = 1.7;
y2 = 0.5;
x3 = 0.9;
y3 = 0.45;

// 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
Mesh.Algorithm = 9;
Mesh.CharacteristicLengthFactor = 3;
Point(1) = {0.1314, 0.0383, 0, cl1};
Point(2) = {0.2745, 0.218, 0, cl1};
Point(3) = {0.3594, 0.3977, 0, cl1};
Point(4) = {0.416, 0.3078, 0, cl1};
Point(5) = {0.4792, 0.2296, 0, cl1};
Point(6) = {0.391, 0.1797, 0, cl1};
Point(7) = {0.2529, 0.1298, 0, cl1};
Point(8) = {0.2978, 0.0715, 0, cl1};
Point(9) = {0.3394, 0.0, 0, cl1};
Point(10) = {0.4759, 0.0549, 0, cl1};
Point(11) = {0.599, 0.1281, 0, cl1};
Point(12) = {0.7321, 0.0566, 0, cl1};
Point(13) = {0.9002, 0.005, 0, cl1};
Point(14) = {0.9451, 0.0765, 0, cl1};
Point(15) = {1.0, 0.1448, 0, cl1};
Point(16) = {0.8453, 0.183, 0, cl1};
Point(17) = {0.7205, 0.2363, 0, cl1};
Point(18) = {0.822, 0.3827, 0, cl1};
Point(19) = {0.8952, 0.5757, 0, cl1};
Point(20) = {0.7887, 0.619, 0, cl1};
Point(21) = {0.7621, 0.612, 0, cl1};
Point(22) = {0.4193, 0.612, 0, cl1};
Point(23) = {0.4276, 0.6606, 0, cl1};
Point(24) = {0.4326, 0.7, 0, cl1};
Point(25) = {0.97, 0.7, 0, cl1};
Point(26) = {0.97, 0.852, 0, cl1};
Point(27) = {0.4493, 0.852, 0, cl1};
Point(28) = {0.4542, 0.9218, 0, cl1};
Point(29) = {0.4576, 0.9983, 0, cl1};
Point(30) = {0.2978, 0.9983, 0, cl1};
Point(31) = {0.2962, 0.9285, 0, cl1};
Point(32) = {0.2928, 0.852, 0, cl1};
Point(33) = {0.0449, 0.852, 0, cl1};
Point(34) = {0.0449, 0.7, 0, cl1};
Point(35) = {0.2762, 0.7, 0, cl1};
Point(36) = {0.1947, 0.4126, 0, cl1};
Point(37) = {0.0, 0.1481, 0, cl1};
Point(38) = {0.0699, 0.0998, 0, cl1};
Point(39) = {0.5923, 0.3278, 0, cl1};
Point(40) = {0.5291, 0.3943, 0, cl1};
Point(41) = {0.4825, 0.464, 0, cl1};
Point(42) = {0.6889, 0.464, 0, cl1};
Point(43) = {0.6439, 0.3894, 0, cl1};

Spline(1) = {1,2,3};
Spline(2) = {3,4,5};
Spline(3) = {5,6,7};
Spline(4) = {7,8,9};
Spline(5) = {9,10,11};
Spline(6) = {11,12,13};
Spline(7) = {13,14,15};
Spline(8) = {15,16,17};
Spline(9) = {17,18,19};
Line(10) = {19,21};
//Line(11) = {20,21};
Line(12) = {21,22};
Spline(13) = {22,23,24};
Line(14) = {24,25};
Line(15) = {25,26};
Line(16) = {26,27};
Spline(17) = {27,28,29};
Line(18) = {29,30};

Spline(19) = {30,31,32};
Line(20) = {32,33};
Line(21) = {33,34};
Line(22) = {34,35};
Spline(23) = {35,36,37};
Spline(24) = {37,38,1};


// furo
Line(25) = {42,41};
Spline(26) = {41,40,39};
Spline(27) = {39,43,42};

Line Loop(1) = {1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24};


Line Loop(2) = {25,26,27};

Plane Surface(1) = {1,2};
Plane Surface(2) = {2};
//Transfinite Line {1} = 3 Using Progression 1;
//Transfinite Line {2} = 2 Using Progression 1;
//Transfinite Line {3} = 2 Using Progression 1;
//Transfinite Line {4} = 2 Using Progression 1;
//Transfinite Line {5} = 3 Using Progression 0.7;
//Transfinite Line {6} = 3 Using Progression 1.42;
//Transfinite Line {7} = 2 Using Progression 1;
//Transfinite Line {8} = 2 Using Progression 1;
//Transfinite Line {9} = 4 Using Progression 0.5;
//Transfinite Line {10} = 2 Using Progression 1;
//Transfinite Line {11} = 0 Using Progression 1;
//Transfinite Line {12} = 3 Using Progression 0.8;
//Transfinite Line {13} = 2 Using Progression 1;
//Transfinite Line {14} = 3 Using Progression 1;
//Transfinite Line {15} = 3 Using Progression 1;
//Transfinite Line {16} = 3 Using Progression 1;
//Transfinite Line {17} = 2 Using Progression 1;
//Transfinite Line {18} = 2 Using Progression 1;
//Transfinite Line {19} = 2 Using Progression 1;
//Transfinite Line {20} = 2 Using Progression 1;
//Transfinite Line {21} = 3 Using Progression 1;
//Transfinite Line {22} = 2 Using Progression 1;
//Transfinite Line {23} = 3 Using Progression 1;
//Transfinite Line {24} = 2 Using Progression 1;
//Transfinite Line {25} = 2 Using Progression 1;
//Transfinite Line {26} = 2 Using Progression 1;
//Transfinite Line {27} = 2 Using Progression 1;


Physical Point(201) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,39,39,40,41,42,43};

//Line Loop(1) = {1};



//Line Loop(2) = {2,3,4,5};
//Line Loop(3) = {6,7,8,9};
//Line Loop(4) = {10,11,12,13};

//Plane Surface(1) = {1,2,3,4};
//Transfinite Line {1} = 16 Using Progression 1;
//Transfinite Line {2,3,4,5,6,7,8,9,10,11,12,13} = 7 Using Progression 1;

//Plane Surface(2) = {2};
//Plane Surface(3) = {3};
//Plane Surface(4) = {4};
//Plane Surface(1) = {1,2,3,4};
//Plane Surface(2) = {2};
//Plane Surface(1) = {1};

//Point{21} In Surface {1};
//Point{22} In Surface {1};
//Point{23} In Surface {1};


//Circle{2} In Surface {1};
//Circle{3} In Surface {1};
//Circle{4} In Surface {1};
//Circle{5} In Surface {1};

Recombine Surface{1};
//Recombine Surface{2};
//Recombine Surface{3};
//Recombine Surface{4};


//Line{2,3,4,5,6,7,8,9,10,11,12,13} In Surface {1};
//Physical Point(201) = {21,22,23};

Physical Line(201) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
Physical Surface(1) = {1};
Physical Surface(2) = {2};