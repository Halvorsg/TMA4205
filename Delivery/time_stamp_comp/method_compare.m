function [k_dep_norm,k_dep_time,k_dep_norm_VC_level,k_dep_time_VC_level,k_dep_norm_VC_pre_smth,k_dep_time_VC_pre_smth,k_dep_norm_VC_post_smth,k_dep_time_VC_post_smth]=method_compare()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
addpath('C:\Users\Markus_F\Documents\GitHub\proj\TMA4205\Delivery');
addpath('C:\Users\Markus_F\Documents\GitHub\proj\TMA4205\GivenCode');
%I0 = double(imread('frame10.png'));
%I1 = double(imread('frame11.png'));


%% setup
acc=10^-8;                  % speciy error improvement
maxit=100;                  % specify a bound on max iterations if no convergence is reached to specified accuracy
loops=0;                    % specify repetition of algorithms to improve evaluation results


%OF vs Vcycle vs PCG
OF_maxit=10000;             %maximum conjugate gradient runthroughs
VC_maxit=1000;              %maximum full V-cycle runthroughs
PCG_maxit=100;              %maximum PCG runthroughs
k_dep_norm=zeros(4,3);
k_dep_time=zeros(4,3);
% V-cycle level compare
maxlevel=5;                 %maximum level depth
k_dep_norm_VC_level=zeros(4,maxlevel);
k_dep_time_VC_level=zeros(4,maxlevel);
% V-cycle smoother compare
max_smth=4;
max_pre_smth=max_smth;      %maximum pre smoothing iterations taken
max_post_smth=max_smth;     %maximum post smoothing iterations taken
smth_step=1;                %increment between smoothing step size
smth_step_dep_vec_size=length([1:smth_step:max_smth]);

k_dep_norm_VC_pre_smth=zeros(4,smth_step_dep_vec_size);
k_dep_time_VC_pre_smth=zeros(4,smth_step_dep_vec_size);
k_dep_norm_VC_post_smth=zeros(4,smth_step_dep_vec_size);
k_dep_time_VC_post_smth=zeros(4,smth_step_dep_vec_size);


%% exemplum
%{
tol=0.01;
lambda=1000;
[I0,I1] = generate_test_images(256);
[norm_PCG_r,time_PCG_stamps] = pic_PCG(I0,I1,tol,PCG_maxit,lambda);
[norm_VC_r,time_VC_stamps]=pic_V_cycle(I0,I1,tol,VC_maxit,5,5,4,lambda);
[norm_OF_r,time_OF_stamps] = pic_OF_cg(I0, I1, tol, OF_maxit,lambda);

plot(time_PCG_stamps,log(norm_PCG_r),time_VC_stamps,log(norm_VC_r),time_OF_stamps,log(norm_OF_r));
%}

%% compute target tolerance
for k=6:9
lambda=2^(k-4);
[I0,I1] = generate_test_images(2^k);
[M,N]=size(I0);
[~,RHS] = rediscretize(I0,I1,M,N,lambda);
tol=norm(RHS)*acc;


%% method compare mode

norm_PCG=0;
time_PCG=0;
norm_VC=0;
time_VC=0;
norm_OF=0;
time_OF=0;


for iterator=0:loops
[norm_PCG_it,time_PCG_stamps_it] = pic_PCG(I0,I1,tol,PCG_maxit,lambda);
norm_PCG=(norm_PCG*iterator+norm_PCG_it)/(iterator+1);
time_PCG=(time_PCG*iterator+time_PCG_stamps_it)/(iterator+1);

[norm_VC_it,time_VC_stamps_it]=pic_V_cycle(I0,I1,tol,VC_maxit,5,5,3,lambda);
norm_VC=(norm_VC*iterator+norm_VC_it)/(iterator+1);
time_VC=(time_VC*iterator+time_VC_stamps_it)/(iterator+1);

[norm_OF_it,time_OF_stamps_it] = pic_OF_cg(I0, I1, tol, OF_maxit,lambda);
norm_OF=(norm_OF*iterator+norm_OF_it)/(iterator+1);
time_OF=(time_OF*iterator+time_OF_stamps_it)/(iterator+1);
end

k_dep_norm(k-5,:)=[norm_OF(end),norm_VC(end),norm_PCG(end)];
k_dep_time(k-5,:)=[time_OF(end),time_VC(end),time_PCG(end)];

% plot(time_PCG,log(norm_PCG),time_VC,log(norm_VC),time_OF,log(norm_OF));


%% V_cycle level compare mode



norm_VC_levels=zeros(maxlevel,1);
time_VC_levels=zeros(maxlevel,1);
norm_VC=0;
time_VC=0;

for level_cnt=1:maxlevel
   for iterator=0:loops
   [norm_VC_it,time_VC_stamps_it]=pic_V_cycle(I0,I1,tol,maxit,5,5,level_cnt,lambda);
    norm_VC=(norm_VC*iterator+norm_VC_it(end))/(iterator+1);
    time_VC=(time_VC*iterator+time_VC_stamps_it(end))/(iterator+1);
   end
   norm_VC_levels(level_cnt)=norm_VC;
   time_VC_levels(level_cnt)=time_VC;
end

k_dep_norm_VC_level(k-5,:)=norm_VC_levels;
k_dep_time_VC_level(k-5,:)=time_VC_levels;

%% V_cycle smoother compare mode

norm_VC_pre_smth=zeros(smth_step_dep_vec_size,1);
time_VC_pre_smth=zeros(smth_step_dep_vec_size,1);
norm_VC_post_smth=zeros(smth_step_dep_vec_size,1);
time_VC_post_smth=zeros(smth_step_dep_vec_size,1);
cnt=0;
for post_smth_cnt=1:smth_step:(max_post_smth*smth_step)
        norm_VC=0;
        time_VC=0;
        for iterator=0:loops
            
            [norm_VC_it,time_VC_stamps_it]=pic_V_cycle(I0,I1,tol,maxit,1,post_smth_cnt,4,lambda);
            norm_VC=(norm_VC*iterator+norm_VC_it(end))/(iterator+1);
            time_VC=(time_VC*iterator+time_VC_stamps_it(end))/(iterator+1);
            
        end
        cnt=cnt+1;
        norm_VC_post_smth(cnt)=norm_VC;
        time_VC_post_smth(cnt)=time_VC;
end
    k_dep_norm_VC_post_smth(k-5,:)=norm_VC_post_smth;
    k_dep_time_VC_post_smth(k-5,:)=time_VC_post_smth;
    
    cnt=0;
for pre_smth_cnt=1:smth_step:(max_pre_smth*smth_step)
        norm_VC=0;
        time_VC=0;
        for iterator=0:loops
            
            [norm_VC_it,time_VC_stamps_it]=pic_V_cycle(I0,I1,tol,maxit,pre_smth_cnt,1,4,lambda);
            norm_VC=(norm_VC*iterator+norm_VC_it(end))/(iterator+1);
            time_VC=(time_VC*iterator+time_VC_stamps_it(end))/(iterator+1);
            
        end
        cnt=cnt+1;
        norm_VC_pre_smth(cnt)=norm_VC;
        time_VC_pre_smth(cnt)=time_VC;
end
    k_dep_norm_VC_pre_smth(k-5,:)=norm_VC_pre_smth;
    k_dep_time_VC_pre_smth(k-5,:)=time_VC_pre_smth;


end
%% plot time against k for OF vs VX vs PCG
f=figure;
plot([6:9],k_dep_time(:,1),[6:9],k_dep_time(:,2),[6:9],k_dep_time(:,3));
g=figure;
plot([1:maxlevel],k_dep_time_VC_level(1,:),[1:maxlevel],k_dep_time_VC_level(2,:),[1:maxlevel],k_dep_time_VC_level(3,:),[1:maxlevel],k_dep_time_VC_level(4,:));
h=figure;
plot([1:smth_step_dep_vec_size],k_dep_time_VC_pre_smth(1,:),[1:smth_step_dep_vec_size],k_dep_time_VC_pre_smth(2,:),[1:smth_step_dep_vec_size],k_dep_time_VC_pre_smth(3,:),[1:smth_step_dep_vec_size],k_dep_time_VC_pre_smth(4,:));
end

