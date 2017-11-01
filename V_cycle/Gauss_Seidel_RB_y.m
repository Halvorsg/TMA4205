function [u,v,norm_r] = Gauss_Seidel_RB_y(u0, v0, lambda, RHS, Syst_mat, M,N, maxit)

%
% [u,v] = Gauss_Seidel_RB(u0, v0, Ix, Iy, rhsu, rhsv, tol, maxit) performs
% the red black method for solving the optical flow problem.
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
% addpath GivenCode
tic
%% Initializing values
% [M,N] = size(Ix); % M rows and N columns
%% Creating B1
B1 = Syst_mat(1:M*N,1:M*N);
B2 = Syst_mat(M*N+1:2*M*N, M*N+1:2*M*N);
IxIy = Syst_mat(1:M*N,M*N+1:2*M*N);
%% Creating IxIy


%% Create full matrix
rhsu = RHS(1:M*N); rhsv = RHS(M*N+1:end);
%Black first, red second
[black,red] = getRedBlack(M,N);
D1 = B1(black,black); D2 = B1(red,red); D3 = B2(black,black); D4 = B2(red,red);
D1_inv = diag(D1).^-1; D2_inv = diag(D2).^-1; D3_inv = diag(D3).^-1; D4_inv = diag(D4).*-1; 
E1 = B1(black,red); F1 = B1(red,black); E2 = B2(black,red); F2 = B2(red,black);
A1 = IxIy(black,black); A2 = IxIy(red,red);
u1 = u0(black); u2 = u0(red); v1 = v0(black); v2 = v0(red); 
b1 = rhsu(black); b2 = rhsu(red); b3 = rhsv(black); b4 = rhsv(red); 

ZEROS = spdiags(zeros(length(D1),1),0,M*N/2,M*N/2);
A = [D1    , E1    , A1    , ZEROS;
     F1    , D2    , ZEROS , A2;
     A1    , ZEROS , D3    , E2;
     ZEROS , A2    , F2    , D4];
cnt = 0;
norm_r = zeros(maxit,1);
b = [b1;b2;b3;b4];
while maxit > cnt
    cnt = cnt+1;
    v2 = D4\(b4 - (F2*v1 + A2*u2));%
    u2 = D2\(b2 - (F1*u1 + A2*v2));%
    v1 = D3\(b3 - (E2*v2 + A1*u1));%
    u1 = D1\(b1 - (E1*u2 + A1*v1));%
    
    norm_r(cnt) = norm(b-A*[u1;u2;v1;v2]);
end
u = [u1;u2];
v = [v1;v2];
u([black,red]) = u;
v([black,red]) = v;

% u = reshape(u,M,N); v = reshape(v,M,N);
%% Visualization
% img = mycomputeColor(u,v); % Have made a small change in this function;
% figure;
% imagesc(img)


end
