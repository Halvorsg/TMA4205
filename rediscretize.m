function [Syst_mat,RHS,int_check] = rediscretize(scaling,I0,I1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Syst_mat=0;
RHS=0;

[M,N] = size(I0);


int_check=isinteger(int16(M)/scaling) & isinteger(int16(N)/scaling);
if int_check

    cut_M=[0:M/scaling-1]*scaling+1;
    cut_N=[0:N/scaling-1]*scaling+1;
    I0=I0(cut_M,cut_N);
    I1=I1(cut_M,cut_N);
    [M,N] = size(I0); 
    
    [It,Ix,Iy] =  dI(I0,I1);
    lambda = 1000;

    %% dxI^2
    Ix = Ix(:);
    e = -lambda*ones(M*N,1);
    diagonal1 = -4*e + Ix.^2;
    B1 = spdiags([e,e,diagonal1,e,e],[-M,-1,0,1,M],M*N,M*N);
    indices1 = (M+1):M:M*N; indices2 = M:M:(M*N-1);
    B1(indices1,indices2) = 0; B1(indices2,indices1) = 0;
    %% dyI^2
    Iy = Iy(:);
    diagonal2 = -4*e + Iy.^2;
    B2 = spdiags([e,e,diagonal2,e,e],[-M,-1,0,1,M],M*N,M*N);
    indices1 = (M+1):M:M*N; indices2 = M:M:(M*N-1);
    B2(indices1,indices2) = 0; B2(indices2,indices1) = 0;
    %% dxI*dyI
    IxIy = spdiags(Ix.*Iy,0,M*N,M*N);
    %% lambda Ah0

    %% dtI*dxI   dtI*dyI
    It = It(:);
    RHS = -[It.*Ix;It.*Iy];


    %% matrix
    Syst_mat = [B1,IxIy;
        IxIy,B2];
    toc
    
end
end

