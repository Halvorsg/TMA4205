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
addpath 'C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205\GivenCode'
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


indices1 = (M+1):M:M*N; indices2 = M:M:(M*N-1);
B1(indices1,indices2) = 0; B1(indices2,indices1) = 0; B2(indices1,indices2) = 0; B2(indices2,indices1) = 0;
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
[x,norm_r,cnt,time] = conjugate_gradient(C,RHS,x0,tol,maxit);
toc
fprintf('Norm of residual = %1.2e \nNumber of steps = %i\n', norm_r(end),cnt)
u = x(1:M*N); v = x(M*N+1:end);
u = reshape(u,M,N); v = reshape(v,M,N);
FLOPS = (1:cnt)*10*(N*M);
%% Visualization
h = figure;
plot(FLOPS,log(norm_r))
ylim([log(10^-8*norm_r(1)),log(norm_r(1))]);
xlim([0,3*10^9])
str = sprintf('Conjugate Gradient with \x03bb = %i',lambda);
title(str)
xlabel('flops')
ylabel('log(||residual||_2)')
saveTightFigure(h,'Convergence_figures/Conjugate Gradient flops')
h = figure;
plot(time,log(norm_r))
str = sprintf('Conjugate Gradient with \x03bb = %i',lambda);
title(str)
xlabel('Time')
ylabel('log(||residual||_2)')
saveTightFigure(h,'Convergence_figures/Conjugate Gradient time')

save('Convergence_figures/time_table_CG','time')
save('Convergence_figures/flops_table_CG','FLOPS')
save('Convergence_figures/res_table_CG','norm_r')
end






