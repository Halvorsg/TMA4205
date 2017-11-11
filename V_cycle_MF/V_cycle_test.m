function V_cycle_test()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205\GivenCode');
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205');

[u0, v0, ~, ~, ~, ~, ~, ~, ~]=initial_values();
I0 = double(imread('frame10.png'));
I1 = double(imread('frame11.png'));
tic
[M,N] = size(I0);
[I0,I1] = imagePreprocessing(I0,I1);
norm1_r=10;
%% Initializing values
step=0;
int_check=isinteger(1);

u1=u0(:);
v1=v0(:);

[Syst_mat,RHS] = rediscretize(I0,I1,M,N,1000);
% tic
A = {};
SMP = {};
cnt = 1;
flag = true;
while norm1_r(end) > 0.01
    tic
    [ u1,v1 ,norm1_r,A,SMP]=subcycle_test(Syst_mat,RHS,1,u1,v1,M,N,1,3,4,norm1_r,A,SMP,flag);                                                 
    plot(log(norm1_r))
    cnt = cnt+1;
    flag = false;
    toc
end
disp(cnt*toc)
u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)


%[u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit] = initial_values();
%[comp_u,comp_v]=OF_cg(u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit);
%norm([u1;v1]-[comp_u;comp_v])
end

