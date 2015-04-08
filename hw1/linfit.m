function [m c] = linfit(points)
% This function calculates the best fit line (y = mx + c) for the given points
% and return m and c
% Assume points are in the format [x1 y1; x2 y2...]
   global mpoints;
   mpoints = points;
   [mc,fval]  = fminunc(@linfit_error, [1 0]);
   m = mc(1);
   c = mc(2);
end

function err = linfit_error(mc)
%error is defined as sum( (yi - (m*xi + c))^2 )
   global mpoints;
   mt = mc(1);
   ct = mc(2); 
   xs = mpoints(:,1);
   ys = mpoints(:,2);
   ys_calc = mt*xs + ct;
   errs = ys - ys_calc;
   err = errs' * errs;
end
