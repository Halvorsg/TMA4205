function method_compare()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
addpath('..\GivenCode');
%I0 = double(imread('frame10.png'));
%I1 = double(imread('frame11.png'));


%% exemplum
%{
[I0,I1] = generate_test_images(640);
[norm_PCG_r,time_PCG_stamps] = pic_PCG(I0,I1,0.01,100);
[norm_VC_r,time_VC_stamps]=pic_V_cycle(I0,I1,0.01,100,5,5,4);
[norm_OF_r,time_OF_stamps] = pic_OF_cg(I0, I1, 0.01, 100);

plot(time_PCG_stamps,log(norm_PCG_r),time_VC_stamps,log(norm_VC_r),time_OF_stamps,log(norm_OF_r));
%}

%% method compare mode
%{
picture_size=512;
loops=10;
tol=0.01;
maxit=100;
norm_PCG=0;
time_PCG=0;
norm_VC=0;
time_VC=0;
norm_OF=0;
time_OF=0;

[I0,I1] = generate_test_images(picture_size);


for iterator=0:loops-1
[norm_PCG_it,time_PCG_stamps_it] = pic_PCG(I0,I1,tol,maxit);
norm_PCG=(norm_PCG*iterator+norm_PCG_it)/(iterator+1);
time_PCG=(time_PCG*iterator+time_PCG_stamps_it)/(iterator+1);

[norm_VC_it,time_VC_stamps_it]=pic_V_cycle(I0,I1,tol,maxit,5,5,4);
norm_VC=(norm_VC*iterator+norm_VC_it)/(iterator+1);
time_VC=(time_VC*iterator+time_VC_stamps_it)/(iterator+1);

[norm_OF_it,time_OF_stamps_it] = pic_OF_cg(I0, I1, tol, maxit);
norm_OF=(norm_OF*iterator+norm_OF_it)/(iterator+1);
time_OF=(time_OF*iterator+time_OF_stamps_it)/(iterator+1);
end

plot(time_PCG,log(norm_PCG),time_VC,log(norm_VC),time_OF,log(norm_OF));
%}

%% V_cycle level compare mode
%{
maxlevel=5;
picture_size=512;
loops=3;
tol=0.01;
maxit=100;

norm_VC_levels=zeros(maxit,maxlevel);
time_VC_levels=zeros(maxit,maxlevel);
norm_VC=0;
time_VC=0;
[I0,I1] = generate_test_images(picture_size);

for level_cnt=1:maxlevel
   for iterator=0:loops
   [norm_VC_it,time_VC_stamps_it]=pic_V_cycle(I0,I1,tol,maxit,5,5,level_cnt);
    norm_VC=(norm_VC*iterator+norm_VC_it)/(iterator+1);
    time_VC=(time_VC*iterator+time_VC_stamps_it)/(iterator+1);
   end
   norm_VC_levels(1:length(norm_VC),level_cnt)=norm_VC;
   time_VC_levels(1:length(time_VC),level_cnt)=time_VC;
end
%}

%% V_cycle smoother compare mode
max_pre_smth=2;
max_post_smth=2;
smth_step=1;
picture_size=512;
loops=2;
tol=0.01;
maxit=100;
[I0,I1] = generate_test_images(picture_size);

norm_VC_smth=zeros(max_pre_smth,max_post_smth);
time_VC_smth=zeros(max_pre_smth,max_post_smth);
norm_VC=0;
time_VC=0;

for post_smth_cnt=1:smth_step:(max_post_smth*smth_step)
    for pre_smth_cnt=1:smth_step:(max_pre_smth*smth_step)
        for iterator=0:loops
            [norm_VC_it,time_VC_stamps_it]=pic_V_cycle(I0,I1,tol,maxit,pre_smth_cnt,post_smth_cnt,4);
            norm_VC=(norm_VC*iterator+norm_VC_it)/(iterator+1);
            time_VC=(time_VC*iterator+time_VC_stamps_it)/(iterator+1);
        end
        norm_VC_smth(post_smth_cnt,pre_smth_cnt)=norm_VC(end);
        time_VC_smth(post_smth_cnt,pre_smth_cnt)=time_VC(end);
    end
end

