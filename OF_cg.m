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
addpath GivenCode
tic
%% Initializing values
[M,N] = size(Ix); % M rows and N columns
%% Creating B1
Ix = Ix(:);
Iy = Iy(:);
e = -lambda*ones(M*N,1);
diagonal1 = -4*e + Ix.^2;
diagonal2 = -4*e + Iy.^2;
B1 = spdiags([e,e,diagonal1,e,e],[-M,-1,0,1,M],M*N,M*N);
B2 = spdiags([e,e,diagonal2,e,e],[-M,-1,0,1,M],M*N,M*N);

% B1hat = spdiags([e,e,diagonal1*0,e,e],[-M,-1,0,1,M],M*N,M*N);
% B2hat = spdiags([e,e,diagonal2*0,e,e],[-M,-1,0,1,M],M*N,M*N);

indices1 = (M+1):M:M*N; indices2 = M:M:(M*N-1);
B1(indices1,indices2) = 0; B1(indices2,indices1) = 0; B2(indices1,indices2) = 0; B2(indices2,indices1) = 0;
% B1hat(indices1,indices2) = 0; B1hat(indices2,indices1) = 0; B2hat(indices1,indices2) = 0; B2hat(indices2,indices1) = 0;
%% Creating IxIy
IxIy = spdiags(Ix.*Iy,0,M*N,M*N);

%% Create full matrix
C = [B1,IxIy;
    IxIy,B2];

%% Creating RHS
RHS = [rhsu(:);rhsv(:)];
toc
%% Solving system Cx = RHS
x0 = [u0(:);v0(:)];
tic
[x,r,cnt,norm_r] = conjugate_gradient(C,RHS,x0,tol,maxit);
toc
fprintf('Norm of residual = %1.2e \nNumber of steps = %i\n', norm(r),cnt)
u = x(1:M*N); v = x(M*N+1:end);
u = reshape(u,M,N); v = reshape(v,M,N);
%% Visualization
img = mycomputeColor(u,v); % Have made a small change in this function;
figure;
imagesc(img)
figure
plot(norm_r)

%%temp
[comp_1,comp_2] = rediscretize(I0,I1,M,N,lambda);
end






