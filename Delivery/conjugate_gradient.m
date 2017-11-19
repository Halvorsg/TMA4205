% The actual itereations in CG
function [x,norm_r,cnt,time] = conjugate_gradient(A,b,x0,tol,maxit)
r = b-A*x0; p = r; x = x0; s = r'*r; s0 = s;
cnt = 1;
norm_r = zeros(maxit,1);
norm_r(cnt) = sqrt(s);
time = zeros(maxit,1);
tic
while sqrt(s)/sqrt(s0) > tol && cnt < maxit
    cnt = cnt+1;
    z = A*p;
    alpha = s/(z'*p);
    x = x + alpha*p;
    r = r-alpha*z;
    s_old = s;
    s = r'*r;
    beta = s/s_old;
    p = r + beta * p;
    norm_r(cnt) = sqrt(s);
    time(cnt) = toc;
end
norm_r = norm_r(1:cnt);
time = time(1:cnt);
end