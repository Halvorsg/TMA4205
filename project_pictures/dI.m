function [dtI,dxI,dyI] = dI(I0,I1)
%% Pre-allocating
[n,m] = size(I0);
dxI0 = zeros(size(I0));
dxI1 = zeros(size(I1));
dyI0 = zeros(size(I0));
dyI1 = zeros(size(I1));
%% dtI
dtI = I1-I0;

%% dxI
for i = 1:m-1
    dxI0(:,i) = I0(:,i+1)-I0(:,i);
    dxI1(:,i) = I1(:,i+1)-I1(:,i);
end
    dxI0(:,m) =  I0(:,m)-I0(:,m-1);
    dxI1(:,m) =  I1(:,m)-I1(:,m-1);
    
    dxI = 1/2*(dxI0+dxI1);
%% dyI
for i = 1:n-1
    dyI0(i,:) = I0(i+1,:)-I0(i,:);
    dyI1(i,:) = I1(i+1,:)-I1(i,:);
end
    dyI0(n,:) =  I0(n,:)-I0(n-1,:);
    dyI1(n,:) =  I1(n,:)-I1(n-1,:);
    
    dyI = 1/2*(dyI0+dyI1);
end