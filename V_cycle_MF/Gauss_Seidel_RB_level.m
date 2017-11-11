function [u,v,A] = Gauss_Seidel_RB_level(u0, v0, RHS, Syst_mat, M,N, maxit,saveMat,level,A, reverse)

%
% [u,v] = Gauss_Seidel_RB(u0, v0, Ix, Iy, rhsu, rhsv, tol, maxit) performs
% the red black method for solving the optical flow problem.
%
% input:
% u0 - initial guess for u
% v0 - initial guess for v
% Ix - x-derivative of the first frame
% Iy - y-derivative of the first frame
% lambda - regularisation parameter
% rhsu - right-hand side in the equation for u
% rhsv - right-hand side in the equation for v
% tol - relative residual tolerance
% maxit - maximum number of iterations
%
% output:
% u - numerical solution for u
% v - numerical solution for v
% norm_r = norm(RHS-Syst_mat*[u0;v0]);
%% Initializing values
switch saveMat
    case true
        B1 = Syst_mat(1:M*N,1:M*N);
        B2 = Syst_mat(M*N+1:2*M*N, M*N+1:2*M*N);
        IxIy = Syst_mat(1:M*N,M*N+1:2*M*N);
        rhsu = RHS(1:M*N); rhsv = RHS(M*N+1:end);
        %Black first, red second
        [black,red] = getRedBlack(M,N);
        D1 = full(diag(B1(black,black))).^-1; D2 = full(diag(B1(red,red))).^-1; 
        D3 = full(diag(B2(black,black))).^-1; D4 = full(diag(B2(red,red))).^-1;
        E1 = B1(black,red); F1 = B1(red,black); E2 = B2(black,red); F2 = B2(red,black);
        A1 = IxIy(black,black); A2 = IxIy(red,red);
        u1 = u0(black); u2 = u0(red); v1 = v0(black); v2 = v0(red); 
        b1 = rhsu(black); b2 = rhsu(red); b3 = rhsv(black); b4 = rhsv(red); 

        A{level} = {D1,D2,D3,D4,E1,E2,F1,F2,A1,A2,black,red,u1,u2,v1,v2};
    case false
        rhsu = RHS(1:M*N); rhsv = RHS(M*N+1:end);
        black = A{level}{11};
        red = A{level}{12};
        switch reverse
            case false
                u1 = A{level}{13};
                u2 = A{level}{14};
                v1 = A{level}{15};
                v2 = A{level}{16};
            case true
                u1 = u0(black); u2 = u0(red); v1 = v0(black); v2 = v0(red); 
        end
        b1 = rhsu(black); b2 = rhsu(red); b3 = rhsv(black); b4 = rhsv(red); 
        D1 = A{level}{1};
        D2 = A{level}{2};
        D3 = A{level}{3};
        D4 = A{level}{4};
        E1 = A{level}{5};
        E2 = A{level}{6};
        F1 = A{level}{7};
        F2 = A{level}{8};
        A1 = A{level}{9};
        A2 = A{level}{10}; 
end

cnt = 0;
switch reverse
    case false
        while maxit > cnt
            cnt = cnt+1;
            u1 = D1.*(b1 - (E1*u2 + A1*v1));%
            u2 = D2.*(b2 - (F1*u1 + A2*v2));%
            v1 = D3.*(b3 - (E2*v2 + A1*u1));%
            v2 = D4.*(b4 - (F2*v1 + A2*u2));%
        end
    case true
        while maxit > cnt
            cnt = cnt+1;
            v2 = D4.*(b4 - (F2*v1 + A2*u2));%
            v1 = D3.*(b3 - (E2*v2 + A1*u1));%
            u2 = D2.*(b2 - (F1*u1 + A2*v2));%
            u1 = D1.*(b1 - (E1*u2 + A1*v1));%
        end
end
% toc
u = [u1;u2];
v = [v1;v2];
u([black;red]) = u;
v([black;red]) = v;
% norm_r = [norm_r,norm([rhsu;rhsv]-Syst_mat*[u;v])];


end
