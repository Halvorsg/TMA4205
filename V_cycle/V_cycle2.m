function V_cycle2()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205\GivenCode');
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205');

[u0, v0, Ix, Iy, ~, rhsu, rhsv, ~, ~] = initial_values();
tic
[M,N] = size(Ix);


%% Initializing values

RHS = [rhsu(:);rhsv(:)];
u1=u0(:);
v1=v0(:);

norm1_r = [];
for i = 1:5
    [u1,v1,norm_r]=subcycle2(0,u1,v1,Ix,Iy,RHS,10,10,1000,4,norm1_r);
    norm1_r = [norm1_r;norm_r];
end

u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)

end

