% The actual itereations in CG
function [x,r,cnt] = conjugate_gradient(A,b,x0,tol,maxit)
r = b-A*x0; p = r; x = x0; s = r'*r;
cnt = 1;
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
end
end