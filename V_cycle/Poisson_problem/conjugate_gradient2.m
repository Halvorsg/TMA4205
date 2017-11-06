% The actual itereations in CG
function [x,r,cnt,norm_r] = conjugate_gradient2(A,b,x0,tol,maxit)
r0 = b-A*x0; r = r0; p = r0; x = x0; s = r0'*r0;
cnt = 1;
r = zeros(length(r0),maxit);
norm_A = [];
r(:,cnt) = r0;
while sqrt(norm(r(:,cnt))) > tol && cnt < maxit
    alpha = dot(r(:,cnt),r(:,cnt))/dot(A*p,p);
    x = x + alpha*p;
    r(:,cnt+1) = r(:,cnt)-alpha*(A*p);
    beta = dot(r(:,cnt+1),r(:,cnt+1))/dot(r(:,cnt),r(:,cnt));
    p = r(:,cnt+1) + beta * p;
    norm_A(cnt) = sqrt(r(:,cnt)'*A*r(:,cnt));
    cnt = cnt+1;

end
figure
norm_r = sqrt(r'*r);
plot(norm_A)
end