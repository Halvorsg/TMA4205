function [ C] = create_cut_mat( M,N  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Cut_mat_row=[1:M*N/2,1:M*N/2,1:M*N/2-M/2,1:M*N/2-M/2];
Cut_mat_column=([1:2:M*N,2:2:M*N,M+1:2:M*N,M+2:2:M*N]);
Cut_mat=sparse(Cut_mat_row,Cut_mat_column,1);
cut_vec=repmat([ones(1,M/2),zeros(1,M/2)],1,N/2);
sp_cut_vec=sparse(1:M*N/2,1:M*N/2,cut_vec);
Cut_mat=sp_cut_vec*Cut_mat;
C=Cut_mat(any(Cut_mat,2),:);
end

