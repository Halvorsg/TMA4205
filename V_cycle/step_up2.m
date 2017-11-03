function [ v1 ] = step_up2( u1,M,N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

u1=reshape(u1,M/2,N/2);
v1=zeros(M,N);
v1(1:2:M,1:2:N)=u1;
v1(2:2:M,1:2:N)=u1;
v1(1:2:M,2:2:N)=u1;
v1(2:2:M,2:2:N)=u1;
v1=v1(:);

v2=zeros(M,N);
v2(1:2:M,1:2:N) = u1;
v2(2:2:M,1:2:N) = 1/2*(u1
end

