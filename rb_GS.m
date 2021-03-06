function [ u1,v1 ] = rb_GS( Syst_mat,RHS,u0,v0,step_cnt,M,N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

u0 = [u0;v0];

black_u0=u0([1:length(u0)/2]*2-1);
red_u0=u0([1:length(u0)/2]*2);

diag_D=diag(Syst_mat);
red_vec=zeros(length(u0)/2,1);
black_vec=zeros(length(u0)/2,1);

for i=1:2M

    if mod(i,2)<0.5
    black_vec(((i-1)*(N/2)+1):(i*N/2))=[1:N/2]*2-1+(i-1)*N;
    red_vec(((i-1)*(N/2)+1):(i*N/2))=[1:N/2]*2+(i-1)*N;
    else
        if mod(i,2)>0.5
            black_vec(((i-1)*(N/2)+1):(i*N/2))=[1:N/2]*2+(i-1)*N;
            red_vec(((i-1)*(N/2)+1):(i*N/2))=[1:N/2]*2-1+(i-1)*N;
        end
    end
    
    
    
end


D=spdiags(diag_D,0,length(diag_D),length(diag_D));
D_inv=spdiags(1./diag_D,0,length(diag_D),length(diag_D));

black_D=D(black_vec,black_vec);
red_D=D(red_vec,red_vec);
black_D_inv=D_inv(black_vec,black_vec);
red_D_inv=D_inv(red_vec,red_vec);

black_RHS=RHS(black_vec);
red_RHS=RHS(red_vec);

U=sparse(triu(Syst_mat));
L=sparse(tril(Syst_mat));
L_U=L+U;

black_L_U=L_U(black_vec,red_vec);
red_L_U=L_U(red_vec,black_vec);

for i=1:step_cnt
red_u0=red_u0+red_D_inv*(red_RHS-(red_L_U)*black_u0+red_D*red_u0);
black_u0=black_u0+black_D_inv*(black_RHS-(black_L_U)*red_u0+black_D*black_u0);
end

u0(black_vec)=black_u0;
u0(red_vec)=red_u0;

u1 = u0(1:M*N); v1 = u0(M*N+1:end);

end

