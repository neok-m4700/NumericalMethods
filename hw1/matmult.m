%This program calculates timings for all 3 versions. It calls the compiled c function
%as well.

function _  = matmult
   for i = 20:20:200
      [t1 t2 t3] = time_it(i);
      printf("%4.0d %0.5g %0.5g %0.5g\n", i, t1, t2, t3);
      fflush(stdout);
   end
end

function [tloop tbuiltin tnative] = time_it(n)
   m1 = generate_matrix(1,n);
   m2 = generate_matrix(2,n);
   tic;
   mult1 = mult_mat_loop(m1,m2);
   tloop = toc;
   tic;
   mult2 = m1*m2;
   tbuiltin = toc;
   [status,cmdout] = system(["./matmult.o " num2str(n)]);
   tnative = str2num(cmdout);
end

function m = generate_matrix(i,n)
%generate a square matrix of nxn with larget element i*i
   a = linspace(0,i,n);
   m = a'*a;
end

function m = mult_mat_loop(m1, m2)
   q = size(m1)(1);
   r = size(m2)(2);
   m = zeros(q,r);
   p = size(m1)(2) ;
   p1 = size(m2)(1);
   assert(p == p1) ;
   for i = 1:q
      for j = 1:r
         for k = 1:p
            m(i,j) += m1(i,k)*m2(k,j);
         end
      end
   end
end


