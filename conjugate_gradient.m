% Conjugate gradient
function [x,r,cnt] = conjugate_gradient(A,b,x0,I0)
r0 = b-x0; r = r0; p = r0; x = x0; s = r0'*r0;
cnt = 1;
while s > eps*(r0'*r0) && cnt < length(A)*2
    z = A*p;
    alpha = s/(z'*p);
    x = x + alpha*p;
    r = r-alpha*z;
    s_old = s;
    s = r'*r;
    beta = s/s_old;
    p = r + beta * p;
    cnt = cnt+1;
end

[M,N] = size(I0);
u = x(1:M*N); v = x(M*N+1:end);
u = reshape(u,size(I0)); v = reshape(v,size(I0));
img = mycomputeColor(u,v); % Have made a small change in this function;
figure;
imagesc(img)
end