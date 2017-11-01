function [ u1,v1 ] = subcycle(level,u0,v0,I0,I1,pre_s,post_s,lambda,max_level)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% getSystem
[M,N] = size(I0);
[Syst_mat,RHS,int_check] = rediscretize(I0,I1,M,N);

int_check=isinteger(int16(M)/(2^level)) && isinteger(int16(N)/(2^level));

%% pre_smooth
[u1,v1]=rb_GS(Syst_mat,RHS,u0,v0,pre_s,M,N);
[u,v] = Gauss_Seidel_RB(u0, v0, Ix, Iy, lambda, rhsu, rhsv, eps, pre_s)

%% initiate recursion
if level<max_level && int_check
[u1]=step_down(u1,M,N);
[v1]=step_down(v1,M,N);
[I0]=step_down_im(I0,M,N);
[I1]=step_down_im(I1,M,N);

[u1,v1]=subcycle(level+1,u1,v1,I0,I1,pre_s,post_s,lambda,max_level);
[u1,v1]=step_up(u1,v1,M,N);
else
    if level==max_level
    %% cg 
    x0=[u1;v1];
    [x,r,cnt]=conjugate_gradient(Syst_mat,RHS,x0,10^-6,500);
    u1 = x(1:M*N); v1 = x(M*N+1:end);
    end
  
end


%% post_smooth
[u1,v1]=rb_GS(Syst_mat,RHS,u1,v1,post_s,M,N);



end

