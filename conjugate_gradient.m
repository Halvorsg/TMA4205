% The actual itereations in CG
function [x,r,cnt,norm_r] = conjugate_gradient(A,b,x0,tol,maxit)
r0 = b-x0; r = r0; p = r0; x = x0; s = r0'*r0;
cnt = 1;
norm_r = [];
while sqrt(s) > tol && cnt < maxit
    z = A*p;
    alpha = s/(z'*p);
    x = x + alpha*p;
    r = r-alpha*z;
    s_old = s;
    s = r'*r;
    beta = s/s_old;
    p = r + beta * p;
    norm_r(cnt) = norm(r);
    cnt = cnt+1;
end

end