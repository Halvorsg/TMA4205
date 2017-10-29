function [ u1,v1 ] = subcycle(level,u0,v0,I0,I1,pre_s,post_s,lambda,max_level)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% getSystem
[Syst_mat,RHS,int_check] = rediscretize(2^level,I0,I1);
[M,N] = size(I0);
M=M/(2^level);
N=N/(2^level);
int_check=isinteger(int16(M)/(2^level)) & isinteger(int16(N)/(2^level));

%% pre_smooth
[u1,v1]=rb_GS(Syst_mat,RHS,u0,v0,pre_s,M,N);

%% initiate recursion
if level<max_level & int_check
[u1,v1]=step_down(u1,v1,M,N);
[u1,v1]=subcycle(level+1,u1,v1,I0,I1,pre_s,post_s,lambda,max_level);
[u1,v1]=step_up(u1,v1,M,N);
end

%% post_smooth
[u1,v1]=rb_GS(Syst_mat,RHS,u1,v1,post_s,M,N);



end

