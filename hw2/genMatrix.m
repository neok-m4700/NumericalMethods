function [ A ] = genMatrix( off_coeff, dims_size)
%GENMATRIX Generate matrix for Heat equation time step 
%  Given a relationship like X(i,j,k) = SumOverL[Y(i+d_li, j+d_lj, k+d_lk)*C_l]
% L is number of terms
% This will return a matrix A, such that AY = X
% Y and X are vectors, that range from 000 to IJK
% IJK = dims_size(1) * dims_size(2) * dims_size(3)
% off_coeff stores the d's and the C's as [ di dj dk C]
% off_coeff has L rows
   IJK = vecSize(dims_size);
   A = zeros(IJK);
   l = length(off_coeff);
   offsets = off_coeff(:, 1:end-1);
   coeffs = off_coeff(:, end);
   for row = 1:IJK
      ijk = idx2ijk(row, dims_size);
      for term = 1:l
         offset = offsets(term, :)';
         ijk_offset = ijk .+ offset;
         col = ijk2idx(ijk_offset, dims_size);
         if all(ijk_offset >= 0)
            A(row,col) = coeffs(term);
         end
      end
   end
end

function [ idx ] = ijk2idx(indices, dims_size)
   idx = 0;
   for d = 1:length(indices)
      idx = idx*dims_size(d) + indices(d);
   end
end

function [ indices ] = idx2ijk( idx, dims_size)
   l = length(dims_size);
   indices = zeros(l,1);
   for d = 1:l
      indices(l-d+1) = mod (idx, dims_size(l-d+1));
      idx = floor( idx / dims_size(l-d+1));
   end
end

function [v] = vecSize(dims_size)
   v = 1;
   for d = 1:length(dims_size)
      v = v*dims_size(d);
   end
end 
