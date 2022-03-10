h = 5;
Mesh.ElementOrder = 1;
Point(1) = {0, 50, 0, h};
Point(2) = {0, 200, 0, h};
Point(3) = {50, 0, 0, h};
Point(4) = {100, 0, 0, h};
Point(5) = {100, 200, 0, h};
Point(6) = {0, 0, 0};
Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};
Circle(5) = {3, 6, 1};
Line Loop(6) = {3, 4, -1, -5, 2};

Plane Surface(7) = {6};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {4};
Physical Surface(4) = {7};
