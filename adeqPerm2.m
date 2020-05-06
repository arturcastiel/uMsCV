

kmap = [1 1 0 0 1; 2 1 0 0 4; 3 4 0 0 1];

ref1 = centelem(:,1) < 0.333;
ref3 = centelem(:,1) > 0.666;
ref2 = ~(ref1 | ref3);

elem(ref1,5) = 1;
elem(ref2,5) = 2;
elem(ref3,5) = 3;
