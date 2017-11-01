function [ v1 ] = step_up( u1,M,N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

u1=reshape(u1,M,N);
v1=zeros(2*M,2*N);
v1(1:2:2*M,1:2:2*N)=u1;
v1(2:2:2*M,1:2:2*N)=u1;
v1(1:2:2*M,2:2:2*N)=u1;
v1(2:2:2*M,2:2:2*N)=u1;

end

