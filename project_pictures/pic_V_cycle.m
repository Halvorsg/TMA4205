function [norm_VC_r,time_stamps]=pic_V_cycle(I0,I1,tol,maxit,pre_s,post_s,max_level)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

time_stamps=zeros(maxit,1);
norm_VC_r=zeros(maxit,1);
tic
[M,N] = size(I0);
u0 = zeros(M,N);
v0 = zeros(M,N);
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
saveMat = true;
time_stamps(1)=toc;
norm_VC_r(1)=norm(RHS);
while norm1_r(end) > tol && cnt<maxit
    tic
    [ u1,v1 ,norm1_r,A,SMP]=subcycle(Syst_mat,RHS,1,u1,v1,M,N,pre_s,post_s,max_level,norm1_r,A,SMP,saveMat);                                                 
    %plot(log(norm1_r))
    cnt = cnt+1;
    saveMat = false;
    time_stamps(cnt) = time_stamps(cnt-1) + toc;
    norm_VC_r(cnt)=norm1_r(end);
%    fprintf('Time passed this round: %f. \t Time passed in total: %f \n',toc,time)
end
time_stamps=time_stamps(1:cnt);
norm_VC_r=norm_VC_r(1:cnt);
%{
u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)

%}
%[u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit] = initial_values();
%[comp_u,comp_v]=OF_cg(u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit);
%norm([u1;v1]-[comp_u;comp_v])
end

