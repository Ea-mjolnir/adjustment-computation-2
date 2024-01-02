%--------------------------------------------------------------------------
%   
%             ADJUSTMENT THEORY II
%          Robust Parameter Estimation  
%              - L1 Adjustment -
% 
%   Author         : Anastasia Pasioti
%   Version        : July 12, 2017
%   Last changes   : July 13, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
%-----------------------Task 2 Final results-------------------------------
%--------------------------------------------------------------------------
%   Observations
%--------------------------------------------------------------------------
% data = 

%Coordinates
% x = 
% y = 

%Vector of observations
% L =

%Number of observations
no_n = length(L);

%Number of unknowns
% no_u = 

%Redundancy
r = no_n-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
% S_LL  

%Theoretical standard deviation
% sigma_0 

%Cofactor matrix of the observations
% Q_LL 

%Weight matrix
% P 

%--------------------------------------------------------------------------
%  L2 - Adjustment
%--------------------------------------------------------------------------
%Design matrix
% A = 

%Normal matrix
% N =

%Vector of right hand side of normal equations
% n = 

%Inversion of normal matrix
% Q_xx = 

%Solution of normal equations
% X_hat = 
      
%Vector of residuals
% v = 

%Vector of adjusted observations
L_hat = L+v;

%Empirical reference standard deviation
s_0 = sqrt(v'*P*v/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknows
s_X = sqrt(diag(S_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));

%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------
%Plot residuals


%Plot results from L2 estimator


%Plot the points





