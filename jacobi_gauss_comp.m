function jacobi_gauss_comp(grid_size)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
[ LHS,  RHS, LHS_diag, LHS_lower, LHS_upper, comp_u] = sparse_create( grid_size );
%ODE_FINITEDIFFERENCE Calculating ODE with basic numerics
%solving -Uxx-2Ux=f with Jacobi and Gausseidel

%% prepare Gaus_Seidel
disc_size = grid_size - 2;
h=1/(grid_size-1);
D=LHS_diag;
E=-LHS_lower;
F=-LHS_upper;

%% Jacobi_matrix
tic;
inv_D=inv(D);
G=inv_D*(E+F);
b=inv_D*RHS;


log_err_new=0;
log_err_old=1;

err=1;
run_count=0;
approx_u=zeros(disc_size,1);
temp_vec=approx_u;


    while abs(log_err_new-log_err_old)>2*eps
        run_count=run_count+1;
        
        temp_vec=G*approx_u+b;
        
        [err,log_err] = calc_err(comp_u,approx_u,temp_vec,disc_size);
        log_err_J(run_count)=log_err;
        
        log_err_old=log_err_new;
        log_err_new=log_err;
        approx_u=temp_vec;
    end
    approx_u_1=approx_u;
    
x_scale_J = 1:length(log_err_J);
time_jacobi=toc;

%% plot both solutions %%%%%%
%x_scale_comp = 1:length(comp_u);
%plot(x_scale_comp,comp_u,x_scale_comp,approx_u);
%end plot

%% Gaus_Seidel 1
%{
tic
D_E_inv=inv(D-E);
G=D_E_inv*F;
b=D_E_inv*RHS;
err=1;
run_count=0;

log_err_new=0;
log_err_old=1;

approx_u=zeros(disc_size,1);
temp_vec=approx_u;

    while abs(log_err_new-log_err_old)>2*eps
        run_count=run_count+1;
        
        temp_vec=G*approx_u+b;
        
        [err,log_err] = calc_err(comp_u,approx_u,temp_vec,disc_size);
        log_err_GS1(run_count)=log_err;
        
        log_err_old=log_err_new;
        log_err_new=log_err;
        
        approx_u=temp_vec;
    end
    
    err=1;
    approx_u_2=approx_u;
 %}    
    tic;
    approx_u=zeros(disc_size,1);
    temp_vec=approx_u;
    run_count=0;
   
      
                temp_vec(1)=(1/LHS(1,1))*(RHS(1)-approx_u(2)*LHS(1,2));
                
                for iterator=2:disc_size-1
                   temp_vec(iterator)=(1/LHS(iterator,iterator)*(RHS(iterator)-temp_vec(iterator-1)*LHS(iterator,iterator-1)-approx_u(iterator+1)*LHS(iterator,iterator+1)));

                end

                temp_vec(disc_size)=(1/LHS(disc_size,disc_size))*(RHS(disc_size)-temp_vec(disc_size-1)*LHS(disc_size,disc_size-1));
    
                
                
    [err,log_err] = calc_err(comp_u,approx_u,temp_vec,disc_size);
        
                approx_u=temp_vec;
        
    err_est=err;
    
    while err>err_est*(10^-6)
        run_count=run_count+1;
        
                temp_vec(1)=(1/LHS(1,1))*(RHS(1)-approx_u(2)*LHS(1,2));
                
                for iterator=2:disc_size-1
                   temp_vec(iterator)=(1/LHS(iterator,iterator)*(RHS(iterator)-temp_vec(iterator-1)*LHS(iterator,iterator-1)-approx_u(iterator+1)*LHS(iterator,iterator+1)));

                end

                temp_vec(disc_size)=(1/LHS(disc_size,disc_size))*(RHS(disc_size)-temp_vec(disc_size-1)*LHS(disc_size,disc_size-1));

                [err,log_err] = calc_err(comp_u,approx_u,temp_vec,disc_size);
                log_err_GS1(run_count)=log_err;
                
               approx_u=temp_vec;
       
            
    end
    
    approx_u_2=approx_u;
    
    plot(1:length(approx_u_2),approx_u_2,1:length(approx_u),approx_u);
    

    
time_Gaus_Seidel_1=toc;

x_scale_GS1 = 1:length(log_err_GS1);

%% Gaus_Seidel 2
%{
tic
D_F_inv=inv(D-F);
G=D_F_inv*E;
b=D_F_inv*RHS;
err=1;
run_count=0;

log_err_new=0;
log_err_old=1;

approx_u=zeros(disc_size,1);
temp_vec=approx_u;

    while abs(log_err_new-log_err_old)>2*eps
        run_count=run_count+1;
        
        temp_vec=G*approx_u+b;
        
        [err,log_err] = calc_err(comp_u,approx_u,temp_vec,disc_size);
        log_err_GS2(run_count)=log_err;
        
        log_err_old=log_err_new;
        log_err_new=log_err;
        
        approx_u=temp_vec;
    end

    approx_u_3=approx_u;
    time_Gaus_Seidel_2=toc;

x_scale_GS2 = 1:length(log_err_GS2);
%}

tic;
    approx_u=zeros(disc_size,1);
    temp_vec=approx_u;
    run_count=0;
   
        
                temp_vec(disc_size)=(1/LHS(disc_size,disc_size))*(RHS(disc_size)-approx_u(disc_size-1)*LHS(disc_size,disc_size-1));
    
                
                for i=1:disc_size-2
                    iterator=disc_size-i;
                   temp_vec(iterator)=(1/LHS(iterator,iterator)*(RHS(iterator)-approx_u(iterator-1)*LHS(iterator,iterator-1)-temp_vec(iterator+1)*LHS(iterator,iterator+1)));

                end
                
                temp_vec(1)=(1/LHS(1,1))*(RHS(1)-temp_vec(2)*LHS(1,2));
                
                
    
    [err,log_err] = calc_err(comp_u,approx_u,temp_vec,disc_size);
        
                approx_u=temp_vec;
    err_est=err;
    
    while err>err_est*(10^-6)
        run_count=run_count+1;
        
                temp_vec(disc_size)=(1/LHS(disc_size,disc_size))*(RHS(disc_size)-approx_u(disc_size-1)*LHS(disc_size,disc_size-1));
    
                
                for i=1:disc_size-2
                    iterator=disc_size-i;
                   temp_vec(iterator)=(1/LHS(iterator,iterator)*(RHS(iterator)-approx_u(iterator-1)*LHS(iterator,iterator-1)-temp_vec(iterator+1)*LHS(iterator,iterator+1)));

                end
                
                temp_vec(1)=(1/LHS(1,1))*(RHS(1)-temp_vec(2)*LHS(1,2));
                
                [err,log_err] = calc_err(comp_u,approx_u,temp_vec,disc_size);
                log_err_GS2(run_count)=log_err;
                
               approx_u=temp_vec;
       
            
    end
    
    approx_u_3=approx_u;
    
    plot(1:length(approx_u_2),approx_u_2,1:length(approx_u),approx_u);
    

    
time_Gaus_Seidel_2=toc;
x_scale_GS2 = 1:length(log_err_GS2);

%% plot results
%plot(x_scale_J,log_err_J,x_scale_GS1,log_err_GS1,x_scale_GS2,log_err_GS2);
plot(disperse_over_intervall(length(comp_u)),comp_u,disperse_over_intervall(length(approx_u_1)),approx_u_1,disperse_over_intervall(length(approx_u_2)),approx_u_2,disperse_over_intervall(length(approx_u_3)),approx_u_3);      

[time_jacobi time_Gaus_Seidel_1 time_Gaus_Seidel_2]


end



