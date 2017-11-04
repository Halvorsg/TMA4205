function [A] = Discretisize(M,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% M rows and N columns
%% Creating B1 and B2
e = -ones(M*N,1);
diagonal1 = -4*e;
A = spdiags([e,e,diagonal1,e,e],[-M,-1,0,1,M],M*N,M*N);
indices1 = (M+1):M:M*N; indices2 = M:M:(M*N-1);
A(indices1,indices2) = 0; A(indices2,indices1) = 0;
%% RHS
end

