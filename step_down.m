function [ u1,v1 ] = step_down(u1,v1,M,N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

u1 = reshape(u1,M,N); v1 = reshape(v1,M,N);
u1 = u1([1:M/2]*2-1,[1:N/2]*2-1);
v1= v1([1:M/2]*2-1,[1:N/2]*2-1);
u1=u1(:);
v1=v1(:);
end

