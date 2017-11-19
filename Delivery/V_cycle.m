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

lambda = 1000;
[Syst_mat,RHS] = rediscretize(I0,I1,M,N,lambda);
norm1_r=norm(RHS);
% tic
A = {};
SMP = {};
maxit = 1000;
cnt = 1;
saveMat = true;
time = zeros(1,maxit);
FLOPS = zeros(1,maxit);
pre_s = 3;
post_s = 3;
max_level = 4;
while norm1_r(end)/norm1_r(1) > 10^-8 && cnt < maxit
    tic
    [ u1,v1 ,~,A,SMP,CG_FLOPS]=subcycle(Syst_mat,RHS,1,u1,v1,M,N,pre_s,post_s,max_level,norm1_r,A,SMP,saveMat,0);                                                 
    cnt = cnt+1;
    saveMat = false;
    norm1_r(cnt) = norm(RHS-Syst_mat*[u1;v1]);
    FLOPS(cnt) = FLOPS(cnt-1) + 7*N*M*(sum((1:max_level).^-1))*(pre_s+post_s)+CG_FLOPS;
    time(cnt) = time(cnt-1) + toc;
    
    %% Time keeping
    fprintf('Time passed this round: %f. \t Time passed in total: %f \n',toc,time(cnt))
    %% Plotting
    %Flops
    figure(1)
    plot(FLOPS(1:cnt),log(norm1_r))
    ylim([log(10^-8*norm1_r(1)),log(norm1_r(1))]);
    xlim([0,3*10^9])
    % Time
    figure(2)
    plot(time(1:cnt),log(norm1_r(1:cnt)))
    ylim([log(10^-8*norm1_r(1)),log(norm1_r(1))]);
    drawnow



end
time = time(1:cnt);
FLOPS = FLOPS(1:cnt);


h = figure(1);
str = sprintf('V-Cycle with \x03bb = %i',lambda);
title(str)
xlabel('flops')
ylabel('log(||residual||_2)')
saveTightFigure(h,'Convergence_figures/V_cycle flops')
h = figure(2);
str = sprintf('V-Cycle with \x03bb = %i',lambda);
title(str)
xlabel('Time')
ylabel('log(||residual||_2)')
saveTightFigure(h,'Convergence_figures/V_cycle time')

save('Convergence_figures/time_table_V_cycle','time')
save('Convergence_figures/flops_table_V_cycle','FLOPS')
save('Convergence_figures/res_table_V_cycle','norm1_r')

u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
img = mycomputeColor(u1,v1); % Have made a small change in this function;
figure;
imagesc(img)


%[u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit] = initial_values();
%[comp_u,comp_v]=OF_cg(u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit);
%norm([u1;v1]-[comp_u;comp_v])
end

