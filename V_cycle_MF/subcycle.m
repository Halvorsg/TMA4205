function [ u1,v1 ,norm1_r] = subcycle(Syst_mat,RHS,level,u0,v0,M,N,pre_s,post_s,max_level,norm1_r)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% getSystem



int_check=isinteger(int16(M)/(2^level)) && isinteger(int16(N)/(2^level));

%% pre_smooth
% [u1,v1]=rb_GS(Syst_mat,RHS,u0,v0,pre_s,M,N);
%norm1_r = [norm1_r;norm(RHS-Syst_mat*[u0;v0])];
%plot(norm1_r,'-*')
[u1,v1,norm_r] = Gauss_Seidel_RB_y(u0, v0, RHS, Syst_mat, M,N, pre_s);
norm1_r = [norm1_r;norm_r];
plot(norm1_r,'-*')
title('Pre-smooth')
% u1 = u0(:);
% v1 = v0(:);

%% initiate recursion
if level<max_level && int_check
[u1]=step_down(u1,M,N);
[v1]=step_down(v1,M,N);
[Syst_mat_2,RHS_2]=step_down_Mat(Syst_mat,RHS,M,N);


[u1,v1,norm1_r]=subcycle(Syst_mat_2,RHS_2,level+1,u1,v1,M/2,N/2,pre_s,post_s,max_level,norm1_r);

[u1]=step_up(u1,M,N);
[v1]=step_up(v1,M,N);
% [I1]=step_up_im(I1,M,N);
else
    if level==max_level
    %% cg 
    x0=[u1;v1];
    [x,r,cnt,norm_r]=conjugate_gradient(Syst_mat,RHS,x0(:),10^-3,100);
    fprintf('Norm of residual = %1.2e \nNumber of steps = %i\n', norm(r),cnt)
    norm1_r = [norm1_r;norm_r'];
    plot(norm1_r,'-*')
    title('After CG')
    u1 = x(1:M*N); v1 = x(M*N+1:end);
    end
  
end


%% post_smooth
[u1,v1,norm_r] = Gauss_Seidel_RB_y(u1, v1, RHS, Syst_mat, M,N, post_s);
r = norm(RHS-Syst_mat*[u1;v1])
norm1_r = [norm1_r;norm_r];
plot(norm1_r);
title('Post smoothing')

end

