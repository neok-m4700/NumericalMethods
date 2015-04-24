function [ A ] = genMatrix2( off_coeffs, dims_size)
%GENMATRIX Generate matrix for Heat equation time step 
%  Given a relationship like X(i,j,k) = SumOverL[Y(i+d_li, j+d_lj, k+d_lk)*C_l]
% L is number of terms
% This will return a matrix A, such that AY = X
% Y and X are vectors, that range from 000 to IJK
% IJK = dims_size(1) * dims_size(2) * dims_size(3)
% off_coeff stores the d's and the C's as [ di dj dk C]
% off_coeff has L rows
nterms = rows(off_coeffs);
N = prod(dims_size); % n_x * n_y * n_z
B = zeros(N, nterms); % columns of B are diagonals of A
d = zeros(2, nterms); % d specifies distance to central diaganoal for each column of B

for t = 1:rows(off_coeffs)
   off_coeff = off_coeffs(t,:);
   coeff = off_coeff(end);
   B(:,t) = coeff;
   offset_nd = off_coeff(1:end-1);
   offset_1d = ijk2idx(offset_nd, dims_size);
   d(t) = offset_1d;
end

A = spdiags(B, d, N, N);
end


function [ idx ] = ijk2idx(indices, dims_size)
   idx = 0;
   for d = 1:length(indices)
      idx = idx*dims_size(d) + indices(d);
   end
end
