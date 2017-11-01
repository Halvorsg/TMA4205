function [u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit] = initial_values();

I0 = double(imread('frame10.png'));
I1 = double(imread('frame11.png'));
% [I0,I1] = generate_test_images(256);
[I0,I1] = imagePreprocessing(I0,I1);
[It,Ix,Iy] =  dI(I0,I1);
lambda = 1000;

rhsu = -Ix.*It;
rhsv = -Iy.*It;

u0 = zeros(size(rhsu));
v0 = zeros(size(rhsv));

tol = 0.001;
maxit = 1000;
end