function [Syst_mat,RHS,int_check] = rediscretize(I0,I1,M,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Syst_mat=0;
RHS=0;

[It,Ix,Iy] =  dI(I0,I1);
% M rows and N columns
%% Creating B1
Ix = Ix(:);
Iy = Iy(:);
e = -lambda*ones(M*N,1);
diagonal1 = -4*e + Ix.^2;
diagonal2 = -4*e + Iy.^2;
B1 = spdiags([e,e,diagonal1,e,e],[-M,-1,0,1,M],M*N,M*N);
B2 = spdiags([e,e,diagonal2,e,e],[-M,-1,0,1,M],M*N,M*N);

B1hat = spdiags([e,e,diagonal1*0,e,e],[-M,-1,0,1,M],M*N,M*N);
B2hat = spdiags([e,e,diagonal2*0,e,e],[-M,-1,0,1,M],M*N,M*N);

indices1 = (M+1):M:M*N; indices2 = M:M:(M*N-1);
B1(indices1,indices2) = 0; B1(indices2,indices1) = 0; B2(indices1,indices2) = 0; B2(indices2,indices1) = 0;
B1hat(indices1,indices2) = 0; B1hat(indices2,indices1) = 0; B2hat(indices1,indices2) = 0; B2hat(indices2,indices1) = 0;
%% Creating IxIy
IxIy = spdiags(Ix.*Iy,0,M*N,M*N);

%% Create full matrix
C = [B1,IxIy;
    IxIy,B2];

end

