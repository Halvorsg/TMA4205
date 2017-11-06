function [ Syst_mat,RHS ] = step_down_Mat(Syst_mat,RHS,M,N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Cut_mat=create_cut_mat(M,N);

B1=Syst_mat(1:M*N,1:M*N);
B2=Syst_mat(M*N+1:end,M*N+1:end);
D1=Syst_mat(1:M*N,M*N+1:end);
D2=Syst_mat(M*N+1:end,1:M*N);

RHS_1=RHS(1:M*N);
RHS_2=RHS(M*N+1:end);

%% clip
B1=Cut_mat*B1*Cut_mat';
B2=Cut_mat*B2*Cut_mat';
D1=Cut_mat*D1*Cut_mat';
D2=Cut_mat*D2*Cut_mat';

Syst_mat=0.25*[B1,D1;D2,B2];
%% columnclip
RHS_1=step_down(RHS_1,M,N);
RHS_2=step_down(RHS_2,M,N);
RHS=[RHS_1;RHS_2];

end

