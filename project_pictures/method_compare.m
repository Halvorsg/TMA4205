function method_compare()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%I0 = double(imread('frame10.png'));
%I1 = double(imread('frame11.png'));

[I0,I1] = generate_test_images(640);

%[norm1_r,time_stamps] = PCG(I0,I1,tol,maxit);
%[norm1_r,time_stamps]=pic_V_cycle(I0,I1,tol,maxit,pre_s,post_s,max_level);
%[norm_r,time_stamps] = OF_cg(I0, I1, tol, maxit);

[norm_PCG_r,time_PCG_stamps] = pic_PCG(I0,I1,0.01,100);
[norm_VC_r,time_VC_stamps]=pic_V_cycle(I0,I1,0.01,100,5,5,4);
[norm_OF_r,time_OF_stamps] = pic_OF_cg(I0, I1, 0.01, 100);
end

