function [m c] = linfit(points)
% This function calculates the best fit line (y = mx + c) for the given points
% and return m and c
% Assume points are in the format [x1 y1; x2 y2...]

   function err = linfit_error(mc)
   %return RMS error. mc = [m c]
      mt = mc(1);
      ct = mc(2); 
      xs = points(:,1);
      ys = points(:,2);
      ys_calc = mt*xs + ct;
      errs = ys - ys_calc;
      err = errs' * errs;

   [mc,fval]  = fminunc(@linfit_error, [1 0])
   m = mc(1);
   c = mc(2);

