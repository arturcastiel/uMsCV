cl1 = 1;
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};
Point(5) = {0.4444444444444444, 0.4444444444444444, 0, cl1};
Point(6) = {0.5555555555555556, 0.4444444444444444, 0, cl1};
Point(7) = {0.5555555555555556, 0.5555555555555556, 0, cl1};
Point(8) = {0.4444444444444444, 0.5555555555555556, 0, cl1};
Point(9) = {0.5, 0.4444444444444444, 0, cl1};
Point(10) = {0.5, 0.5555555555555556, 0, cl1};
Point(11) = {0.5, 0, 0, cl1};
Point(12) = {0.5, 1, 0, cl1};
Line(1) = {1, 11};
Line(2) = {11, 2};
Line(3) = {2, 3};
Line(4) = {3, 12};
Line(5) = {12, 4};
Line(6) = {4, 1};
Line(7) = {5, 9};
Line(8) = {9, 6};
Line(9) = {6, 7};
Line(10) = {7, 10};
Line(11) = {10, 8};
Line(12) = {8, 5};
Line(13) = {9, 11};
Line(14) = {10, 12};
Line Loop(16) = {1, -13, -7, -12, -11, 14, 5, 6};
Plane Surface(16) = {16};
Line Loop(18) = {2, 3, 4, -14, -10, -9, -8, 13};
Plane Surface(18) = {18};
Physical Point(19) = {1, 2, 3, 4, 11, 12};
Physical Point(20) = {5, 6, 7, 8, 9, 10};
Physical Line(21) = {1, 2, 3, 4, 5, 6};
Physical Line(22) = {7, 8, 9, 10, 11, 12};
Physical Line(23) = {13, 14};
Physical Surface(24) = {16, 18};
