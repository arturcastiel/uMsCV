Mesh.Algorithm = 2;
Mesh.Algorithm =7;
//Mesh.CharacteristicLengthFactor = 1;
cl__1 = 1.125;
Point(1) = {0, 0, 0, 1.125};
Point(2) = {0.5, 0, 0, 1.125};

Point(3) = {1, 0, 0, 1.125};
Point(4) = {1, 0.5, 0, 1.125};
Point(5) = {1, 1, 0, 1.125};

Point(6) = {0.5, 1, 0, 1.125};

Point(7) = {0, 1, 0, 1.125};
Point(8) = {0, 0.5, 0, 1.125};



Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {2, 8};
Line(10) = {3, 7};
Line(11) = {4, 6};


//Line Loop(1) = {1, 2, 3, 4,5,6,7,8};
Line Loop(1) = {1, 9, 8};
Line Loop(2) = {2, 10, 7,-9};
Line Loop(3) = {3,11,6,-10};
Line Loop(4) = {4,5,-11};

//Transfinite Line {10} = 18 Using Progression 1;

//Transfinite Line {1,4,5, 8} = 20 Using Progression 1;
//Transfinite Line {2, 3,6,7} = 20 Using Progression 1;
//Transfinite Line {9,11} = 23 Using Progression 1;
//Transfinite Line {10} = 36 Using Progression 1;

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};


Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};

Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};

Physical Point(201) = {1, 2, 3, 4};
Physical Line(201) = {1, 2, 3, 4,5,6,7,8};