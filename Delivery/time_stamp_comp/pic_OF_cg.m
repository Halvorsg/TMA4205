function [norm_r,time_stamps] = pic_OF_cg(I0, I1, tol, maxit,lambda)
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
tic
% [I0,I1] = generate_test_images(256);
[I0,I1] = imagePreprocessing(I0,I1);
[It,Ix,Iy] =  dI(I0,I1);
rhsu = -Ix.*It;
rhsv = -Iy.*It;
u0 = zeros(size(rhsu));
v0 = zeros(size(rhsv));
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
setup_time=toc;
%% Solving system Cx = RHS
x0 = [u0(:);v0(:)];
time_stamps=zeros(maxit,1);
time_stamps(1)=setup_time;

[x,r,cnt,norm_r,time_stamps] = pic_conjugate_gradient(C,RHS,x0,tol,maxit,time_stamps);

%fprintf('Norm of residual = %1.2e \nNumber of steps = %i\n', norm(r),cnt)
%u = x(1:M*N); v = x(M*N+1:end);
%u = reshape(u,M,N); v = reshape(v,M,N);
%% Visualization
%{
img = mycomputeColor(u,v); % Have made a small change in this function;
figure;
imagesc(img)
figure
plot(norm_r)

%%temp
%[comp_1,comp_2] = rediscretize(I0,I1,M,N,lambda);
%}
end






