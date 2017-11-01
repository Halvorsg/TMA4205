function [ v1 ] = step_down(u1,M,N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

u1 = reshape(u1,M,N);
v1 = 1/4*(u1(1:2:M,1:2:N)+u1(2:2:M,1:2:N)+u1(1:2:M,2:2:N)+u1(2:2:M,2:2:N));
%u1=u1(:);

end

