function [ Syst_Mat,RHS ] = step_up_Mat( Syst_Mat,RHS,M,N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

B1=Syst_mat(1:M*N,1:M*N);
B2=Syst_mat(M*N+1:end,M*N+1:end);
D1=Syst_mat(1:M*N,M*N+1:end);
D2=Syst_mat(M*N+1:end,1:M*N);

RHS_1=RHS(1:M*N);
RHS_2=RHS(M*N+1:end);


end

