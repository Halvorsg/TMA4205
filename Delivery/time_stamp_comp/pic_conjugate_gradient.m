% The actual itereations in CG
function [x,r,cnt,norm_r,time_stamps] = conjugate_gradient_2(A,b,x0,tol,maxit,time_stamps)
r = b-A*x0; p = r; x = x0; s = r'*r;
cnt = 1;
norm_r=zeros(maxit,1);
norm_r(1)=norm(b);

while sqrt(s) > tol && cnt < maxit
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
time_stamps=time_stamps(1:cnt);
end