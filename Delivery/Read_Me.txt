functions supposed to run as mains

function V_cycle()
runs a V_Cycle to accuracy improvemen 10^-8 for the presented reallife images

function [x,r,cnt] = PCG(tol,maxit)
runs the PCG method to given accuracy with maximal iterations max_it and gives back the solution,
the residual and the amount of iterations taken

function [u,v] = OF_cg(u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit)
runs the CG method for the reallife images. starting values can be obtained with [...]=initial_values();


function[...]=method_compare();
gives back arrayys which conclude runtime of certain algorithms tested on images from
the test image function with sze 2^k, k=6:9
criteria for the tests is specified in the function.
main criteria: 	acc 	- accuracy of computation
		maxit	- maximum amount of iterations per method
		loops	- repetition per method 

functions used (own functions)

function [norm_r,time_stamps] = pic_OF_cg(I0, I1, tol, maxit,lambda)
function [norm_PCG_r,time_stamps] = pic_PCG(I0,I1,tol,maxit,lambda)
function [norm_VC_r,time_stamps]=pic_V_cycle(I0,I1,tol,maxit,pre_s,post_s,max_level,lambda)
track norm and computation time for CG, PCG and VC method. 

function [ u1,v1 ,norm1_r,A,SMP] = pic_subcycle(Syst_mat,RHS,level,u0,v0,M,N,pre_s,post_s,max_level,norm1_r,A,SMP,flag)
contains the actual V-cylce for pic_V_cycle

function [x,norm_r,cnt] = conjugate_gradient(A,b,x0,tol,maxit)
performs the CG method on Ax=b


function [C] = create_cut_mat(M,N)
creates the restriction operator

function [dtI,dxI,dyI] = dI(I0,I1)
computes the descretization of the pictures

function [u,v,A] = Gauss_Seidel_RB_level(u0, v0, RHS, Syst_mat, M,N, maxit,saveMat,level,A,reverse)
performs red_black_Gauss Seidel computation with respect to the secial problem

function [black,red] = getRedBlack(M,N)
gets the coloring for the points with respect to our special problem

function [I0,I1] = imagePreprocessing(I0,I1)
applies a gaussian filter to the images


function [u0, v0, Ix, Iy, lambda, rhsu, rhsv, tol, maxit] = initial_values();
gets the inital values to run the CG main function


function [Syst_mat,RHS] = rediscretize(I0,I1,M,N,lambda)
computes the rystem matrix and rigth hand side for specified puctures

function [ v1 ] = step_down(u1,M,N)
restriction of a vector to coarser level

function [ Syst_mat,RHS ] = step_down_Mat_level(Syst_mat,RHS,M,N,saveMat)
restricton of Matrix and right hand side to coarser level

function [ v1 ] = step_up( u1,M,N )
prolongation of vector

function [ u1,v1 ,norm1_r,A,SMP,FLOPS] = subcycle(Syst_mat,RHS,level,u0,v0,M,N,pre_s,post_s,max_level,norm1_r,A,SMP,flag,FLOPS)
actual V-cycle recursion


functions used (external)

function [I1,I2] = generate_test_images(n) 

function img = mycolorwheel(n)

function img = mycomputeColor(u,v)


