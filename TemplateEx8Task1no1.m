%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION
%        Variance component estimation  
% 
%   Author         : Anastasia Pasioti
%   Version        : June 21, 2017
%   Last changes   : June 20, 2018
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format long g;
%---------------------------Task 1.1---------------------------------------
%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Load all files
dist = load('Distances.txt');
dir = load('Directions.txt');

%Vector of observations
L = [dist(:,3); dir(:,3)*pi/200];    %Convert to [rad]

%Gauss-Krueger coordinates for fixed points [m]
y6 = 5317651.428;
x6 = 4968940.373;
y9 = 5324162.853;
x9 = 4970922.160;

%Initial values for new points [m]
y1 = 5314698.13;
x1 = 4965804.18;
y15 = 5320448.85;
x15 = 4962997.53;

%Initial values for orientation unknowns
w1 = 0;
w6 = 0;
w9 = 0;
w15 = 0;

%Initial values for unknowns
X_0 = [y1 x1 y15 x15 w1 w6 w9 w15]';

%Number of observations
no_n = length(L);

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_dist = 0.1;             %[m]
s_dir = 0.001*pi/200;     %Convert to [rad]

s_LL = [s_dist^2*ones(length(dist),1); s_dir^2*ones(length(dir),1)];
S_LL = diag(s_LL);

%Theoretical reference standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);
P_dist = P(1:5,1:5); %change: weight for dist
P_dir = P(6:14,6:14); %change: weigh for dir    


%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-5;
delta = 10^-9;
max_x_hat = 10^Inf;

%Number of iterations
iteration = 0;

%Initialization A
A = zeros(no_n,no_u);

%Iteration

while max_x_hat>epsilon || Check2>delta

    %Vector of distances
    
    L_0(1)=dis(y6,x6,y1,x1);
    L_0(2)=dis(y9,x9,y1,x1);
    L_0(3)=dis(y9,x9,y6,x6);
    L_0(4)=dis(y15,x15,y1,x1);
    L_0(5)=dis(y15,x15,y9,x9);
    
    %Vector of directions
    
    L_0(6)=direction(y1,x1,y6,x6,w1);
    L_0(7)=direction(y1,x1,y15,x15,w1);
    L_0(8)=direction(y6,x6,y1,x1,w6);
    L_0(9)=direction(y6,x6,y9,x9,w6);
    L_0(10)=direction(y9,x9,y15,x15,w9);
    L_0(11)=direction(y9,x9,y1,x1,w9);
    L_0(12)=direction(y9,x9,y6,x6,w9);
    L_0(13)=direction(y15,x15,y1,x1,w15);
    L_0(14)=direction(y15,x15,y9,x9,w15);
    
    %Vector of reduced observations
    
    l = L-L_0';

    %Design matrix
    A(1,1)=ds_dy_to(y6,x6,y1,x1);
    A(1,2)=ds_dx_to(y6,x6,y1,x1);
    
    A(2,1)=ds_dy_to(y9,x9,y1,x1);
    A(2,2)=ds_dx_to(y9,x9,y1,x1);
    
    A(4,1)=ds_dy_to(y15,x15,y1,x1);
    A(4,2)=ds_dx_to(y15,x15,y1,x1);
    A(4,3)=ds_dy_from(y15,x15,y1,x1);
    A(4,4)=ds_dx_from(y15,x15,y1,x1);

    A(5,3)=ds_dy_from(y15,x15,y9,x9);
    A(5,4)=ds_dx_from(y15,x15,y9,x9);
    
    A(6,1)=dr_dy_from(y1,x1,y6,x6);
    A(6,2)=dr_dx_from(y1,x1,y6,x6);
    A(6,5)=-1;
    
    A(7,1)=dr_dy_from(y1,x1,y15,x15);
    A(7,2)=dr_dx_from(y1,x1,y15,x15);
    A(7,3)=dr_dy_to(y1,x1,y15,x15);
    A(7,4)=dr_dx_to(y1,x1,y15,x15);
    A(7,5)=-1;
    
    A(8,1)=dr_dy_to(y6,x6,y1,x1);
    A(8,2)=dr_dx_to(y6,x6,y1,x1);
    A(8,6)=-1;
    
    A(9,6)=-1;
    
    A(10,3)=dr_dy_to(y9,x9,y15,x15);
    A(10,4)=dr_dx_to(y9,x9,y15,x15);
    A(10,7)=-1;
    
    A(11,1)=dr_dy_to(y9,x9,y1,x1);
    A(11,2)=dr_dx_to(y9,x9,y1,x1);
    A(11,7)=-1;
    
    A(12,7)=-1;
    
    A(13,1)=dr_dy_to(y15,x15,y1,x1);
    A(13,2)=dr_dx_to(y15,x15,y1,x1);
    A(13,3)=dr_dy_from(y15,x15,y1,x1);
    A(13,4)=dr_dx_from(y15,x15,y1,x1);
    A(13,8)=-1;
    
    A(14,3)=dr_dy_from(y15,x15,y9,x9);
    A(14,4)=dr_dx_from(y15,x15,y9,x9);
    A(14,8)=-1;
    
    %Normal matrix
    N = A'*P*A;

    %Vector of right hand side of normal equations
    n = A'*P*l;

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx = inv(N);
    
    %Solution of the normal equations
    x_hat = Q_xx*n;
    
    %Adjusted unknowns
    X_hat = X_0+x_hat;
    
    %Update
    X_0 = X_hat;
    
    y1 = X_0(1);
    x1 = X_0(2);
    y15 = X_0(3);
    x15 = X_0(4);
    w1 = X_0(5);
    w6 = X_0(6);
    w9 = X_0(7);
    w15 = X_0(8);
    
    %Check 1
    max_x_hat = max(abs(x_hat));
	
	%Vector of residuals
    v = A*x_hat-l;

    %Objective function
    vTPv = v'*P*v;

    %Vector of adjusted observations
    L_hat = L+v;

    % distances
    
    Phi_X_hat(1) = dis(y6,x6,y1,x1);
    Phi_X_hat(2) = dis(y9,x9,y1,x1);
    Phi_X_hat(3) = dis(y9,x9,y6,x6);
    Phi_X_hat(4) = dis(y15,x15,y1,x1);
    Phi_X_hat(5) = dis(y15,x15,y9,x9);
    
    % directions
    
    Phi_X_hat(6) = direction(y1,x1,y6,x6,w1);   
    Phi_X_hat(7) = direction(y1,x1,y15,x15,w1);
    Phi_X_hat(8) = direction(y6,x6,y1,x1,w6);
    Phi_X_hat(9) = direction(y6,x6,y9,x9,w6);
    Phi_X_hat(10) = direction(y9,x9,y15,x15,w9);
    Phi_X_hat(11) = direction(y9,x9,y1,x1,w9);
    Phi_X_hat(12) = direction(y9,x9,y6,x6,w9);
    Phi_X_hat(13) = direction(y15,x15,y1,x1,w15);
    Phi_X_hat(14) = direction(y15,x15,y9,x9,w15);

    %Check 2
    Check2 = max(abs(L_hat-Phi_X_hat'));
    
    %Update number of iterations
    iteration = iteration+1;

end


%Sub-vectors of residuals
v_dist = v(1:5,1);          %sub-vector of v
v_dir = v(6:14,1);          %sub-vector of v
v_dir_gon = v_dir*200/pi;    %Convert to [gon]

%Adjusted observations - Convert to [gon]
L_hat_gon = L_hat(6:14)*200/pi;   

% Final Check
if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);    %a posteriori 

%VC matrix of the adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));
s_X_gon = s_X(5:8,1)*200/pi;        %Convert to [gon]

%Cofactor matrix of the adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of the adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));
s_L_hat_gon = s_L_hat(6:14)*200/pi;  %Convert to [gon]

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%Compute Q_vv*P
Q_vvP = Q_vv*P;
QvvP_dist = Q_vvP(1:5,1:5); 
QvvP_dir = Q_vvP(6:14,6:14); 


%Compute trace for sub-matrices Q_vvP_dist and Q_vvP_dir
r_dist = trace(QvvP_dist);
r_dir = trace(QvvP_dir);

%A posteriori for each group of observations
s_o_dist = sqrt((v_dist'*P_dist*v_dist)/r_dist);
s_o_dir = sqrt((v_dir'*P_dir*v_dir)/r_dir);

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));
s_v_gon = s_v(6:14,1)*200/pi;      %Convert to [gon]

%Results for the unknowns
y1 = X_0(1);
x1 = X_0(2);
y15 = X_0(3);
x15 = X_0(4);

X_final = [y1 x1 y15 x15 w1 w6 w9 w15]';

%Convert to [gon] and check the quadrants
gon = X_final(5:8,1)*200/pi;
gon(1) = gon(1)+400;
gon(2) = gon(2)+400;
gon(4) = gon(4)+400;








