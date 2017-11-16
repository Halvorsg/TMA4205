function [x,r,cnt] = PCG(tol,maxit)
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205\GivenCode');
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205\V_cycle_MF');
addpath('C:\Users\halvo\Documents\MATLAB\TMA4205\TMA4205');
%% Initializing values
I0 = double((imread('frame10.png'))); % Frame 1
I1 = double((imread('frame11.png'))); % Frame 2
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
[ zu1,zv1 ,~,A,SMP]=subcycle_test(Syst_mat,RHS,1,u1(:),v1(:),M,N,3,3,4,norm1_r,A,SMP,flag);
flag = false;
z = [zu1;zv1];
 p = z; s = z'*r;
ZEROS = zeros(size(u1));
cnt = 1;
while sqrt(norm(s)) > tol && cnt < maxit
    tic
    cnt = cnt+1;
    Ap = Syst_mat*p;
    alpha = s/(Ap'*p);
    x = x + alpha*p;
    r = r-alpha*Ap;
    [ zu1,zv1 ,~,A,SMP]=subcycle_test(Syst_mat,r,1, ZEROS, ZEROS ,M,N,3,3,4,norm1_r,A,SMP,flag);
    z = [zu1;zv1];
    s_old = s;
    s = r'*z;
    beta = s/s_old;
    p = z + beta * p;
    norm1_r = [norm1_r;norm(r)];
    toc
end
disp(cnt*toc)
r = RHS-Syst_mat*x;
plot(log(norm1_r))
u1 = x(1:length(u1)); v1 = x(length(v1)+1:end);
u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)
end