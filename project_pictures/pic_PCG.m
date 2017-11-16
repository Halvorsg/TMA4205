function [norm_PCG_r,time_stamps] = pic_PCG(I0,I1,tol,maxit)
addpath('..\GivenCode');
%% Initializing values
time_stamps=zeros(maxit,1);
norm_PCG_r=zeros(maxit,1);
[M,N] = size(I0);                   % Size of the frames
[I0,I1] = imagePreprocessing(I0,I1);% Smooting the pictures
u1 = zeros(size(I0(:))); v1 = zeros(size(I0(:)));   % Initial guess for solution
[Syst_mat,RHS] = rediscretize(I0,I1,M,N,1000);      % Discretized system
A = {};                             % Saving values for Gauss-Seidel method
SMP = {};                           % Saving the matrices for each level
flag = true;                        % Indicator weather to save values or get values
%% PCG
x = [u1;v1]; r = RHS-Syst_mat*x;
norm1_r= norm(r);
[ zu1,zv1 ,~,A,SMP]=subcycle(Syst_mat,RHS,1,u1(:),v1(:),M,N,3,3,4,norm1_r,A,SMP,flag);
flag = false;
z = [zu1;zv1];
 p = z; s = z'*r;
ZEROS = zeros(size(u1));
cnt = 1;
time_stamps(1)=toc;
norm_PCG_r(1)=norm(RHS);
while sqrt(norm(s)) > tol && cnt < maxit
    tic
    cnt = cnt+1;
    Ap = Syst_mat*p;
    alpha = s/(Ap'*p);
    x = x + alpha*p;
    r = r-alpha*Ap;
    [ zu1,zv1 ,~,A,SMP]=subcycle(Syst_mat,r,1, ZEROS, ZEROS ,M,N,3,3,4,norm1_r,A,SMP,flag);
    z = [zu1;zv1];
    s_old = s;
    s = r'*z;
    beta = s/s_old;
    p = z + beta * p;
    norm1_r = [norm1_r;norm(r)];
    time_stamps(cnt) = time_stamps(cnt-1) + toc;
    norm_PCG_r(cnt)=norm1_r(end);
%    fprintf('Time passed this round: %f. \t Time passed in total: %f \n',toc,time)
end
r = RHS-Syst_mat*x;
norm1_r = [norm1_r;norm(r)];
time_stamps=time_stamps(1:cnt);
norm_PCG_r=norm_PCG_r(1:cnt);
%{
plot(log(norm1_r))
u1 = x(1:length(u1)); v1 = x(length(v1)+1:end);
u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)
%}
end