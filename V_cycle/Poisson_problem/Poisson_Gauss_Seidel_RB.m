function [u,norm_r] = Poisson_Gauss_Seidel_RB(u0, RHS, Syst_mat, M,N, maxit)
%% Creating B1
B1 = Syst_mat(1:M*N,1:M*N);
%% Creating IxIy


%% Create full matrix
rhsu = RHS(1:M*N);
%Black first, red second
[black,red] = getRedBlack(M,N);
D1 = B1(black,black); D2 = B1(red,red);
E1 = B1(black,red); F1 = B1(red,black);

u1 = u0(black); u2 = u0(red);
b1 = rhsu(black); b2 = rhsu(red);

A = [D1    , E1;
     F1    , D2]; 

cnt = 0;
norm_r = zeros(maxit,1);
b = [b1;b2];
while maxit > cnt
    cnt = cnt+1;
    u1 = D1\(b1 - E1*u2);%
    u2 = D2\(b2 - F1*u1);%
    
    norm_r(cnt) = norm(b-A*[u1;u2]);
end
u = [u1;u2];
u([black,red]) = u;

%% Visualization
% img = mycomputeColor(u,v); % Have made a small change in this function;
% figure;
% imagesc(img)


end
