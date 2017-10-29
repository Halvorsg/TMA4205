function [u,v] = OF_cg(u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit)
%
% [u,v] = OF_cg(u0, v0, Ix, Iy, rhsu, rhsv, tol, maxit) performs
% the CG method for solving the optical flow problem.
%
% input:
% u0 - initial guess for u
% v0 - initial guess for v
% Ix - x-derivative of the first frame
% Iy - y-derivative of the first frame
% lambda - regularisation parameter
% rhsu - right-hand side in the equation for u
% rhsv - right-hand side in the equation for v
% tol - relative residual tolerance
% maxit - maximum number of iterations
%
% output:
% u - numerical solution for u
% v - numerical solution for v

[M,N] = size(Ix);
%% Creating tri-diagonal part of matrix
Ix = Ix(:);
Iy = Iy(:);
e = -lambda*ones(M*N,1);
diagonal1 = -4*e + Ix.^2;
diagonal2 = -4*e + Iy.^2;
B1 = spdiags([e,e,diagonal1,e,e],[-M,-1,0,1,M],M*N,M*N);
B2 = spdiags([e,e,diagonal2,e,e],[-M,-1,0,1,M],M*N,M*N);
indices1 = (M+1):M:M*N; indices2 = M:M:(M*N-1);
B1(indices1,indices2) = 0; B1(indices2,indices1) = 0; B2(indices1,indices2) = 0; B2(indices2,indices1) = 0;

%% Creating off-diagonal part of matrix
IxIy = spdiags(Ix.*Iy,0,M*N,M*N);
%% Creating final matrix
C = [B1,IxIy;
    IxIy,B2];
%% Creating RHS
RHS = [rhsu(:);rhsv(:)];
%% Solving system
x0 = [u0(:);v0(:)];
[X,~,~] = conjugate_gradient(C,RHS,x0,Ix);

u = reshape(X(1:M*N),size(Ix)); v = reshape(X(M*N+1:2*M*N),size(Ix));







