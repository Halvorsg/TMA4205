function [x,r,cnt] = PCG(tol,maxit)
addpath('..\GivenCode');
%% Initializing values
I0 = double((imread('frame10.png'))); % Frame 1
I1 = double((imread('frame11.png'))); % Frame 2
[M,N] = size(I0);                   % Size of the frames
[I0,I1] = imagePreprocessing(I0,I1);% Smooting the pictures
u1 = zeros(size(I0(:))); v1 = zeros(size(I0(:)));   % Initial guess for solution
lambda = 1000;
[Syst_mat,RHS] = rediscretize(I0,I1,M,N,lambda);      % Discretized system
A = {};                             % Saving values for Gauss-Seidel method
SMP = {};                           % Saving the matrices for each level
flag = true;                        % Indicator weather to save values or get values
%% PCG
x = [u1;v1]; r = RHS-Syst_mat*x;
[ zu1,zv1 ,~,A,SMP,~]=subcycle(Syst_mat,RHS,1,u1(:),v1(:),M,N,3,3,4,norm(r),A,SMP,flag,0);
flag = false;
z = [zu1;zv1];
 p = z; s = z'*r;
ZEROS = zeros(size(u1));
cnt = 1;
time = zeros(maxit,1);
FLOPS = zeros(maxit,1);
norm1_r = zeros(maxit,1);
norm1_r(1) = norm(r);
max_level = 4;
pre_smooth = 1;
post_smooth = 1;
while norm1_r(cnt)/norm1_r(1) > tol && cnt < maxit
    tic
    cnt = cnt+1;
    Ap = Syst_mat*p;
    alpha = s/(Ap'*p);
    x = x + alpha*p;
    r = r-alpha*Ap;
    [ zu1,zv1 ,~,A,SMP,CG_FLOPS]=subcycle(Syst_mat,r,1, ZEROS, ZEROS ,M,N,pre_smooth,post_smooth,max_level,norm1_r,A,SMP,flag,0);
    z = [zu1;zv1];
    s_old = s;
    s = r'*z;
    beta = s/s_old;
    p = z + beta * p;
    FLOPS(cnt) = FLOPS(cnt-1) + 7*N*M*(sum((1:max_level).^-1))*(pre_smooth+post_smooth)+CG_FLOPS + 11*N*M;
    norm1_r(cnt) = norm(r);
    time(cnt) = time(cnt-1) + toc;
    fprintf('Time passed this round: %f. \t Time passed in total: %f \n',toc,time(cnt))
end
%% Saving figures and vectors
h = figure(1);
plot(FLOPS(1:cnt),log(norm1_r(1:cnt)))
ylim([log(10^-8*norm1_r(1)),log(norm1_r(1))]);
xlim([0,3*10^9])
str = sprintf('Preconditioned Conjugate Gradient with \x03bb = %i',lambda);
title(str)
xlabel('flops')
ylabel('log(||residual||_2)')
saveTightFigure(h,'Convergence_figures/PCG flops')

h = figure(2);
plot(time(1:cnt),log(norm1_r(1:cnt)))
ylim([log(10^-8*norm1_r(1)),log(norm1_r(1))]);
str = sprintf('Preconditioned Conjugate Gradient with \x03bb = %i',lambda);
title(str)
xlabel('Time')
ylabel('log(||residual||_2)')
saveTightFigure(h,'Convergence_figures/PCG time')

save('Convergence_figures/time_table_PCG','time')
save('Convergence_figures/flops_table_PCG','FLOPS')
save('Convergence_figures/res_table_PCG','norm1_r')

%% Visualizing the solution
u1 = x(1:length(u1)); v1 = x(length(v1)+1:end);
u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
h = figure;
imagesc(img)
axis off
saveTightFigure(h,'Convergence_figures/result_lambda_1000')
end