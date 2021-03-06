function C = setUpSystem(I0,I1)%u0,v0,Ix,Iy,lambda,rhsu,rhsv,tol,maxit)
addpath GivenCode
tic
%% Initializing values
[I0,I1] = imagePreprocessing(I0,I1);
[It,Ix,Iy] =  dI(I0,I1);
lambda = 1000;
[M,N] = size(Ix);
%% Creating B1
Ix = Ix(:);
e = -lambda*ones(M*N,1);
diagonal1 = -4*e + Ix.^2;
B1 = spdiags([e,e,diagonal1,e,e],[-M,-1,0,1,M],M*N,M*N);
indices1 = (M+1):M:M*N; indices2 = M:M:(M*N-1);
B1(indices1,indices2) = 0; B1(indices2,indices1) = 0;
%% Creating B2
Iy = Iy(:);
diagonal2 = -4*e + Iy.^2;
B2 = spdiags([e,e,diagonal2,e,e],[-M,-1,0,1,M],M*N,M*N);
indices1 = (M+1):M:M*N; indices2 = M:M:(M*N-1);
B2(indices1,indices2) = 0; B2(indices2,indices1) = 0;
%% Creating IxIy
IxIy = spdiags(Ix.*Iy,0,M*N,M*N);

%% Creating [Ix*It;Iy*It] (RHS for u and v)
It = It(:);
RHS = -[It.*Ix;It.*Iy];

%% Create full matrix
C = [B1,IxIy;
    IxIy,B2];
toc
%% Solving system
tic
[x,r,cnt] = conjugate_gradient(C,RHS,zeros(size(RHS)),I0);
toc
% tic
% X = C\RHS;
% toc
% u = X(1:M*N); v = X(M*N+1:end);
% u = reshape(u,size(I0)); v = reshape(v,size(I0));
% img = mycomputeColor(u,v); % Have made a small change in this function;
% figure;
% imagesc(img)
end
