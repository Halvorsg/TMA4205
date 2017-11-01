function V_cycle()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('..\GivenCode');

[u0, v0, ~, ~, ~, ~, ~, ~, ~]=initial_values();
I0 = double(imread('frame10.png'));
I1 = double(imread('frame11.png'));
tic
[M,N] = size(I0);
[I0,I1] = imagePreprocessing(I0,I1);
% 
% [It,Ix,Iy] =  dI(I0,I1);

%% Initializing values
[I0,I1] = imagePreprocessing(I0,I1);
step=0;
int_check=isinteger(1);

u1=u0(:);
v1=v0(:);


[u1,v1]=subcycle(1,u1,v1,I0,I1,10,10,1000,1);
               %(level,u0,v0,I0,I1,pre_s,post_s,lambda,max_level)


u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)

end

