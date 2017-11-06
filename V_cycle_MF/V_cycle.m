function V_cycle()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('C:\Users\Markus_F\Documents\GitHub\proj\TMA4205\GivenCode');
addpath('C:\Users\Markus_F\Documents\GitHub\proj\TMA4205');

[u0, v0, ~, ~, ~, ~, ~, ~, ~]=initial_values();
I0 = double(imread('frame10.png'));
I1 = double(imread('frame11.png'));
tic
[M,N] = size(I0);
[I0,I1] = imagePreprocessing(I0,I1);
norm1_r=0;
%% Initializing values
step=0;
int_check=isinteger(1);

u1=u0(:);
v1=v0(:);

[Syst_mat,RHS] = rediscretize(I0,I1,M,N,1000);

norm1_r = [];
for i = 1:10
    [ u1,v1 ,norm1_r]=subcycle(Syst_mat,RHS,1,u1,v1,M,N,20,5,4,norm1_r);
                                                 
    %norm1_r = [norm1_r;norm_r];
end

u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)

end

