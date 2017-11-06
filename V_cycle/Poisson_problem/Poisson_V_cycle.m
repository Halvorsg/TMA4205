function Poisson_V_cycle(M,N)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205\GivenCode');
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205\V_cycle');
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205');

tic
u0 = zeros(M*N,1);


%% Initializing values

RHS = ones(M*N,1);
u1=u0(:);

norm1_r = [];
for i = 1:5
    [ u1 ,norm1_r] = Poisson_subcycle(0,u0,RHS,10,10,norm1_r,M,N);
    norm1_r = [norm1_r;norm_r];
end

u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)

end

