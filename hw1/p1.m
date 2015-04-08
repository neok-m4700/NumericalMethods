% Let square = |  8  x1  x2 |
%              | x3   5  x4 |
%              | x5   1  x6 |
%
% System of equations:
% Let sum of each row, column and diagonal be y
% First row: 8 + x1 + x2 = y
% rewrite as: x1 + x2 - y = -8
% Similarly for all rows and columns and diagonals
% Total sum = x1 + x2 + .. x6 + 8 + 5 + 1 = 45
% rewrite as: x1 + .. + x6 = 31

% We get the following matrix:

A = [ 1 1 0 0 0 0 -1 -8;
      0 0 1 1 0 0 -1 -5;
      0 0 0 0 1 1 -1 -1;
      0 0 1 0 1 0 -1 -8;
      1 0 0 0 0 0 -1 -6;
      0 1 0 1 0 1 -1  0;
      0 0 0 0 0 1 -1 -13;
      0 1 0 0 1 0 -1 -5;
      1 1 1 1 1 1  0 31; ];
      
B = rref(A)

% Solving this we get:
%  |  8   9  -2 |
%  | -5   5  15 |
%  | 12   1   2 |
% This meets all the constraints. 



