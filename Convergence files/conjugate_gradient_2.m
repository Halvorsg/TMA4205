% The actual itereations in CG
function [x,r,cnt,norm_r,time_stamps] = conjugate_gradient_2(A,b,x0,tol,maxit,time_stamps)
r0 = b-A*x0; r = r0; p = r; x = x0; s = r'*r;
cnt = 1;
norm_r=zeros(maxit,1);
norm_r(1)=norm(b);
time_stamps = zeros(maxit,1);
while sqrt(s)/norm(r0) > tol && cnt < maxit
    tic
    z = A*p;
    alpha = s/(z'*p);
    x = x + alpha*p;
    r = r-alpha*z;
    s_old = s;
    s = r'*r;
    beta = s/s_old;
    p = r + beta * p;
    
    cnt = cnt+1;
    time_stamps(cnt)=time_stamps(cnt-1)+toc;
    norm_r(cnt)=sqrt(s);
end
norm_r=norm_r(1:cnt);
FLOPS = (1:length(norm_r))*10*length(r);
time_stamps=time_stamps(1:cnt);
end