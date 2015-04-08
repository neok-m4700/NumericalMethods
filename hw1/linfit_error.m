function err = linfit_error(points)

% Return the Sum of square of errors
% Error for any points is (yi - (m*xi + c))


xs = points(:,1);
ys = points(:,2);

ys_calc = m*xs + c;
errs = ys - ys_calc;
errsq = err.*err;
err = sum(errsq);

