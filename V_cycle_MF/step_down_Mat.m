function [ Syst_mat,RHS ] = step_down_Mat(Syst_mat,RHS,M,N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Cut_mat_row=[1:M*N/2,1:M*N/2,1:M*N/2-M/2,1:M*N/2-M/2];
Cut_mat_column=([1:2:M*N,2:2:M*N,M+1:2:M*N,M+2:2:M*N]);
Cut_mat=sparse(Cut_mat_row,Cut_mat_column,1);
cut_vec=repmat([ones(1,M),zeros(1,M)],1,N/4);
Cut_mat=Cut_mat.*transpose(cut_vec);
Cut_mat( ~any(Cut_mat,2), : ) = [];


B1=Syst_mat(1:M*N,1:M*N);
B2=Syst_mat(M*N+1:end,M*N+1:end);
D1=Syst_mat(1:M*N,M*N+1:end);
D2=Syst_mat(M*N+1:end,1:M*N);

RHS_1=RHS(1:M*N);
RHS_2=RHS(M*N+1:end);

%% rowclip


%% columnclip


end

