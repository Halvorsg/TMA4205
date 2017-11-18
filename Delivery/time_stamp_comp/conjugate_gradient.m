% The actual itereations in CG
function [x,norm_r,cnt] = conjugate_gradient(A,b,x0,tol,maxit)
r = b-A*x0; p = r; x = x0; s = r'*r; s0 = s;
cnt = 1;
norm_r = zeros(maxit,1);
norm_r(cnt) = sqrt(s);
while sqrt(s) > tol && cnt < maxit
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
end
norm_r = norm_r(1:cnt);
end