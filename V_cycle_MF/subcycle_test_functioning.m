function [ u1,v1 ,norm1_r,A] = subcycle_test_functioning(Syst_mat,RHS,level,u0,v0,M,N,pre_s,post_s,max_level,norm1_r,A,flag)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% getSystem



int_check=isinteger(int16(M)/(2^level)) && isinteger(int16(N)/(2^level));



%% initiate recursion
if level<max_level && int_check

%% pre_smooth
% [u1,v1]=rb_GS(Syst_mat,RHS,u0,v0,pre_s,M,N);
%norm1_r = [norm1_r;norm(RHS-Syst_mat*[u0;v0])];
%plot(norm1_r,'-*')
[u1,v1,A] = Gauss_Seidel_RB_level(u0, v0, RHS, Syst_mat, M,N, pre_s,flag,level,A);
residual_u1_v1=RHS-Syst_mat*[u1;v1];
norm1_r = [norm1_r;norm(residual_u1_v1)];
% plot(log(norm1_r),'-*')
% title('Pre-smooth')
% u1 = u0(:);
% v1 = v0(:);
    
    
%[u0]=step_down(u1,M,N);
%[v0]=step_down(v1,M,N);
% u0=zeros(length(u0)/2,1);
% v0=zeros(length(v0)/2,1);
% tic
[Syst_mat_2,RHS_2]=step_down_Mat(Syst_mat,residual_u1_v1,M,N);
u0=zeros(length(RHS_2)/2,1);
v0=zeros(length(RHS_2)/2,1);
% fprintf('Time for step_down_Mat %f \n',toc)


[up_u,up_v,norm1_r,A]=subcycle_test(Syst_mat_2,RHS_2,level+1,u0,v0,M/2,N/2,pre_s,post_s,max_level,norm1_r,A,flag);

[up_u]=step_up(up_u,M,N);
[up_v]=step_up(up_v,M,N);

u1=u1+up_u;
v1=v1+up_v;
% [I1]=step_up_im(I1,M,N);

%% post_smooth
[u1,v1,A] = Gauss_Seidel_RB_level(u1, v1, RHS, Syst_mat, M,N, pre_s,flag,level,A);
r = norm(RHS-Syst_mat*[u1;v1]);
norm1_r = [norm1_r;norm(r)];
% plot(log(norm1_r));
% title('Post smoothing')


else
    if level==max_level
    %% cg 
%     x0=[u0;v0];
    %[x,r,cnt,norm_r]=conjugate_gradient(Syst_mat,RHS,x0(:),10^-3,100);
    x=Syst_mat\RHS;
    %fprintf('Norm of residual = %1.2e \nNumber of steps = %i\n', norm(r),cnt)
    %norm1_r = [norm1_r;norm_r'];
    %plot(norm1_r,'-*')
    %title('After CG')
    u1 = x(1:M*N); v1 = x(M*N+1:end);
    end
  
end



end

