M = 10;
e = ones(M,1);
A = spdiags([e,2*e,e],-1:1,M,M)
x = 1:2:M;
y = 2:2:M;
A = A([x,y],[x,y]);
full(A)




M = 11;N = 10;
e = ones(M*N,1);
diagonal1 = -4*e;
diagonal2 = -4*e;
A = spdiags([e,e,diagonal1,e,e],[-M,-1,0,1,M],M*N,M*N);
indices1 = [M+1:M:M*N]; indices2 = [M:M:M*N-1];
A(indices1,indices2) = 0; A(indices2,indices1) = 0;
if mod(M,2) == 1
    black = 1:2:M*N;
    red = 2:2:M*N;
    spy(A([black,red],[black,red]))
else
X = reshape(1:M*N,M,N)';
tic %Fastest
red = zeros(ceil(M*N/2),1);
redCNT = 0;
for i = 1:2:M
    red(redCNT+1:redCNT + length(X(2:2:end,i))) = X(2:2:end,i);
    redCNT = redCNT + length(X(2:2:end,i));
end
for i = 2:2:M
    red(redCNT+1:redCNT + length(X(1:2:end,i))) = X(1:2:end,i);
    redCNT = redCNT + length(X(1:2:end,i));
end
toc 
red = sort(red);
black = ([1:2:N*M]+[2:2:M*N])'-red;

spy(A([black,red],[black,red]))

end 

b = rand(length(A),1);
b1 = b(black);
b2 = b(red);
u1 = zeros(size(b1));
u2 = zeros(size(b2));
D1 = A(black,black); D2 = A(red,red);
E = A(black,red) ; F = A(red,black);
cnt = 0;
maxiter = 1000;
r = 1;
norm_r = zeros(maxiter,1);
while r > 0.001 && cnt < maxiter
    cnt = cnt+1;
    u1 = D1\(b1-(E*u2));
    u2 = D2\(b2-(F*u1));
    norm_r(cnt) = norm([b1;b2] - A([black,red],[black,red])*[u1;u2]);
    r = norm_r(cnt);
end
re = length(red); bl = length(black);
plot(log(norm_r))
B = A([black,red],[black,red]);
c = zeros(size(b));
X = zeros(size(A));
X([black,red],[black,red]) = B([1:bl,re+1:end],[1:bl,re+1:end]);
disp(isequal(A,X))

   
    
    