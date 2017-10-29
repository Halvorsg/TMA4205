function [ u1,v1 ] = step_up( u1_p,v1_p,M,N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
u1_p = reshape(u1_p,M/2,N/2);
v1_p = reshape(v1_p,M/2,N/2);
u1 = zeros(2*M,2*N);
v1 = zeros(2*M+2,2*N+2);
u1([1:M/2]*2,[1:N/2]*2)=u1_p;
v1([1:M/2]*2,[1:N/2]*2)=v1_p;
v1=conv2(u1,1/8*[2 4 2;4 8 4;2 4 2]);
v1=conv2(v1,1/8*[2 4 2;4 8 4;2 4 2]);
u1=u1(2:M+1,2:N+1);
v1=v1(2:M+1,2:N+1);
u1=u1(:);
v1=v1(:);

end

