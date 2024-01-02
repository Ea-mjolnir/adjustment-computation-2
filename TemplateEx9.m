%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION
%        Gauss-Helmert model - Part I  
% 
%   Author         : Anastasia Pasioti
%   Version        : June 27, 2017
%   Last changes   : June 27, 2017
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Observations...THE MEASUREMENTS apo thn ekfwnhsh!!
x = [1.0 2.0 3.0 4.0]';     %[m]
y = [0.1 1.1 1.8 2.4]';     %[m]

%Vector of observations
l = [x; y];

%Initial values for the unknowns
a = 1; %slope of the line.....taking two points from the measurements and solving the 2*2 system we find a=1
b = 0; %y-intersection

%Vector of Unknowns
X_0 = [a b]';

%Number of observations/points
no_n = length(l)/2; %Not all the observations!!! THE NUMBER OF POINTS - NUMBER OF CONDITION EQUATIONS!!!
                    % just 4 point
%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_xy = ([2 1 4 2])*0.01;     %[m]

S_ll = diag([s_xy.^2 s_xy.^2]); %std is given but for the S_ll we want the VARIANCES..that's why s_xy^2

%Theoretical reference standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_ll = 1/sigma_0^2*S_ll;

%Weight matrix
P = inv(Q_ll);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-8; %SMALL VALUE
delta = 10^-12; % the break off condition for linearization...our problem NON-LINEAR!!
max_psi = 10^Inf; %MAXIMUM NUMBER OF Ø:CONDITION EQUATION

%Number of iterations
iteration = 0;

%Initialization A and B
A = zeros(no_n,no_u); %...with respect to the unknown parameters.
B = zeros(no_n,no_n*2); %...with respect to the residuals.

%Initial values for the residuals
vx = zeros(no_n,1); %we give the value zero as vector!!!
vy = zeros(no_n,1);

while max_psi>epsilon || Check2>delta
    
    %Psi function..CONDITION EQUATIONS....HAVE TO BE UPDATED ALL TIME
     psi = (y + vy)-a*(x + vx) -b;
    
   %Matrices with the elements from the Jacobian matrices J1, J2
     A = [-x-vx -ones(no_n,1)];

     B = [-a*eye(no_n) eye(no_n)];
    
    %Vector of misclosures
     w = -B*[vx; vy] + psi -0;  %C=0
    
    %Normal matrix
     N_ext = [B*Q_ll*B'  A; 
                 A' zeros(2,2)]; %zeros(2,2) because we have two unknows!!
     
    %Vector of right hand side of normal equations
     n_ext = [-w; 0; 0]; %two 0 because we have two unknows!!
    
    %Inversion of normal matrix / Cofactor matrix of the unknowns..............................
     Q_xx_ext = inv(N_ext);
    
    %Solution of the normal equations...........................................................X_hat
     x_hat = Q_xx_ext*n_ext;
    
    %Vector of residuals...until here the residuals have the value 0!!
     v_new = Q_ll*B'*x_hat(1:no_n);
    
    %Update
    v = v_new;
    % 4 residuals for x and 4 residuals for y
    vx = v(1:end/2);
    vy = v(end/2+1:end); %otherwise v(5:8);
    
    %Update for the Unknowns
    X_0 = X_0 + x_hat(end-no_u+1:end); %We take the last 2 values of vector [k X]...IN GAUSS HELMERT MODEL THE X is below that Lagrange Multi
                                       % last 2 valuse because 2 unkown
                                       % must be updated
    a = X_0(1);
    b = X_0(2);
    
    %Check 1
    max_psi = max(abs(psi));
	
	%vTPv function...aposteriori std of the Adjustement
	vTPv = v'*P*v;

	%Vector of adjusted observations
	L_hat = l+v; %residuals have been updated v=v_new
	
	%Check 2
	Check2 = max(abs(L_hat(5:8,1)-a*L_hat(1:4,1)-b)); % I create again the condition equations by the adjusted observations!!!
	                                                  % see the obs.eqation
                                                      % 5to8 is y, 1to4 is x
	%Update number of iterations
    iteration = iteration+1;
    
end

% Final Check
if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

Q_22 = Q_xx_ext(end-no_u+1:end, end-no_u+1:end);
Q_11 = Q_xx_ext(1:no_n,1:no_n);

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);    %a posteriori 

%VC matrix of adjusted unknowns
C_xx = -s_0^2*Q_xx_ext(end-no_u+1:end, end-no_u+1:end);

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(C_xx));

%Cofactor matrix of the residuals..................................
Q_vv = Q_ll*B'*Q_xx_ext(1:no_n,1:no_n)*B*Q_ll;

%VC matrix of residuals
C_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(C_vv));

%Cofactor matrix of adjusted observations................................
Q_ll_hat = Q_ll-Q_vv;

%VC matrix of adjusted observations
S_ll_hat = s_0^2*Q_ll_hat;

%Standard deviation of the adjusted observations
s_l_hat = sqrt(diag(S_ll_hat));

%Results for the unknowns...parameters of the line!!
a = X_0(1,1);
b = X_0(2,1); 












