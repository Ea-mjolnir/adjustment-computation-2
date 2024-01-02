%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION
%          Robust Parameter Estimation  
%              - L1 Adjustment -
% 
%   Author         : Anastasia Pasioti
%   Version        : July 12, 2017
%   Last changes   : July 12, 2017
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
%-----------------------Task 1---------------------------------------------
%--------------------------------------------------------------------------
%   Observations 
%--------------------------------------------------------------------------
%Load data
data = load('testseries.txt');

%Vector of observations
L = data  %we don't have function model here...not spesific observations!! NOT SPECIFIC PROBLEM!

%Number of observations
no_n = length(L);

%Number of unknowns 
no_u = 1;    %becuase there is nothing else , just pick 1 

%Redundancy
r = no_n-no_u;

%--------------------------------------------------------------------------
%  Initial stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
S_LL = eye(no_n); 

%Theoretical standard deviation
sigma_0 = 1;

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  L1 - Adjustment
%--------------------------------------------------------------------------
%break-off conditions
 epsilon = 1e-14; %As initial values we use very small values---close to 0. If sth goes wrong we increase this values
 c =eps;
 max_v_hat = 10e10;
 
%Initialization
 vk1 = ones(no_n,1); %we need to initialize the residual vector
 
%Number of iterations
iteration = 0;
 
while max_v_hat>epsilon

    %Update of the residuals
      vk = vk1;

    %Design matrix
    A = ones(no_n,1);

    %Normal matrix
     N = A'*P*A;

    %Vector of right hand side of normal equations
     n = A'*P*L;

    %Inversion of normal matrix
     Q_xx = inv(N); 

    %Solution of normal equations
     x_hat = Q_xx*n;

    %Calculation of the new residuals
     vk1 = A*x_hat-L;

    %Update of the weight matrix
    %c is the break off condition that we gave
     P = diag(1./(abs(vk1)+c)); %STEP 3 of the task 1 - is given in the algorithm
   
    %Check
     max_v_hat = max(abs(vk1-vk)); %Algorirthm step 4
   
    %Update number of iterations
    iteration = iteration+1;

end

v = vk1
P_final = diag(P) %the final weight matrix...of the last iteration










