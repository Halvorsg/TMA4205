function V_cycle()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit]=initial_values();
I0 = double(imread('frame10.png'));
I1 = double(imread('frame11.png'));
addpath GivenCode
tic
[M,N] = size(I0);

[It,Ix,Iy] =  dI(I0,I1);

%% Initializing values
[I0,I1] = imagePreprocessing(I0,I1);
step=0;
int_check=isinteger(1);

u0=u0(:);
v0=v0(:);

[u1,v1]=subcycle(0,u0,v0,I0,I1,10,10,1000,4);
               %(level,u0,v0,I0,I1,pre_s,post_s,lambda,max_level)


u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
end

