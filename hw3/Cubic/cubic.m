function _ = cubic(n)
   clf
   hold on
   points = gen_points(n);
   B = [ones(n-2,1) ones(n-2, 1)*4 ones(n-2,1)];
   A = spdiags(B, [-1 0 1], n-2, n-2);
   xs = points(:,1);
   ys = points(:,2);
   h = points(2,1) - points(1,1);
   b = (ys(1:n-2) - 2*ys(2:n-1) + ys(3:n)) * 6 / (h*h);
   M = A\b;
   M = [0 ; M ; 0];
   as = (M(2:end) - M(1:end-1))/(6*h);
   bs = M(1:end-1)/2;
   cs = (ys(2:end) - ys(1:end-1))/h - (M(2:end) + 2*M(1:end-1))*(h/6);
   ds = ys(1:end-1);
   
   res = 100;
   pxs = linspace(xs(1:end-1), xs(2:end), res);
   pys = zeros(n-1,res);
   for i=1:n-1
      for j=1:res
         pys(i,j) = as(i)*(pxs(i,j) - xs(i))^3 + bs(i)*(pxs(i,j) - xs(i))^2 + cs(i)*(pxs(i,j) - xs(i)) + ds(i);
      end
   end
   plot(xs, ys, "bo", 'MarkerSize', 5)
   plot(pxs, pys, "r.")
   hold off
end

function points = gen_points(n)
   x = 1:n;
   y = rand(1,n);
   points = [x' y'];
end
   
