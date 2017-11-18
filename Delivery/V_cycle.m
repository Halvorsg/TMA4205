function V_cycle()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('..\GivenCode');

[u0, v0, ~, ~, ~, ~, ~, ~, ~]=initial_values();
 I0 = double(imread('frame10.png'));
 I1 = double(imread('frame11.png'));

%[I0,I1] = generate_test_images(640);
tic
[M,N] = size(I0);
[I0,I1] = imagePreprocessing(I0,I1);

%% Initializing values
step=0;
int_check=isinteger(1);

u1=u0(:);
v1=v0(:);

[Syst_mat,RHS] = rediscretize(I0,I1,M,N,1000);
norm1_r=norm(RHS);
% tic
A = {};
SMP = {};
maxit = 200;
cnt = 1;
saveMat = true;
time = 0;
FLOPS = zeros(1,maxit);
while norm1_r(end)/norm1_r(1) > 10^-10 && cnt < maxit
    tic
    [ u1,v1 ,~,A,SMP,CG_FLOPS]=subcycle(Syst_mat,RHS,1,u1,v1,M,N,3,3,4,norm1_r,A,SMP,saveMat,0);                                                 
    cnt = cnt+1;
    FLOPS(cnt) = FLOPS(cnt-1) + 5*N*M*(sum((1:4).^-1))*(3*3)+CG_FLOPS;
    norm1_r(cnt) = norm(RHS-Syst_mat*[u1;v1]);
    figure(1)
    plot(FLOPS(1:cnt),log(norm1_r))
    ylim([log(10^-10*norm1_r(1)),log(norm1_r(1))]);
    xlim([0,3*10^9])
    drawnow
    saveMat = false;
    time = time + toc;
    fprintf('Time passed this round: %f. \t Time passed in total: %f \n',toc,time)
end
title('log(||r||) as a function of flops')
u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)


%[u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit] = initial_values();
%[comp_u,comp_v]=OF_cg(u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit);
%norm([u1;v1]-[comp_u;comp_v])
end

