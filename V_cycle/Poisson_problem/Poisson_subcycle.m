function [ u1 ,norm1_r] = Poisson_subcycle(level,u0,RHS,pre_s,post_s,norm1_r,M,N)
max_level = 4;
%% getSystem
Syst_mat = Discretisize(M,N);

int_check=isinteger(int16(M)/(2^level)) && isinteger(int16(N)/(2^level));

%% pre_smooth
% [u1,v1]=rb_GS(Syst_mat,RHS,u0,v0,pre_s,M,N);
norm1_r = [norm1_r;norm(RHS-Syst_mat*u0)];
plot(norm1_r,'-*')
[u1,norm_r] = Poisson_Gauss_Seidel_RB(u0, RHS, Syst_mat, M,N, pre_s);
norm1_r = [norm1_r;norm_r];
plot(norm1_r,'-*')
title('Pre-smooth')
% u1 = u0(:);
% v1 = v0(:);

%% initiate recursion
if level<max_level && int_check
[u1]=step_down(u1,M,N);
[rhsu]=step_down(RHS,M,N);
RHS_new = [rhsu];

    [ u1 ,norm1_r] = Poisson_subcycle(level+1,u1,RHS_new,10,10,norm1_r,M/2,N/2);

[u1]=step_up(u1,M,N);

% [I0]=step_up_im(I0,M,N);
% [I1]=step_up_im(I1,M,N);
else
    if level==max_level
    %% cg 
    x0=[u1];
    [x,r,cnt,norm_r]=conjugate_gradient(Syst_mat,RHS,x0(:),10^-3,1000);
    fprintf('Norm of residual = %1.2e \nNumber of steps = %i\n', norm(r),cnt)
    norm1_r = [norm1_r;norm_r'];
%     plot(norm1_r,'-*')
    title('After CG')
    u1 = x;
    end
  
end


%% post_smooth
[u1,norm_r] = Poisson_Gauss_Seidel_RB(u0, RHS, Syst_mat, M,N, pre_s);
r = norm(RHS-Syst_mat*[u1])
norm1_r = [norm1_r;norm_r];
plot(norm1_r);
title('Post smoothing')

end

