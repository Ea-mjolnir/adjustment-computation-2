%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION
%    Quality Assessment of Adjustment Results  
% 
%   Author         : Anastasia Pasioti
%   Version        : June 14, 2017
%   Last changes   : June 11, 2018
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format long g;
%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
%Load all files
dist = load('Distances.txt');
dir = load('Directions.txt');

%Vector of observations
L = [dist(:,3); dir(:,3)*pi/200];    %Convert to [rad]

%Gauss-Krueger coordinates for new points [m]
y1 = 5314698.13;
x1 = 4965804.18;
y6 = 5317651.428;
x6 = 4968940.373;
y9 = 5324162.853;
x9 = 4970922.160;
y15 = 5320448.85;
x15 = 4962997.53;

%Initial values for orientation unknowns
w1 = 0;
w6 = 0;
w9 = 0;
w15 = 0;

%Initial values for unknowns
X_0 = [y1 x1 y6 x6 y9 x9 y15 x15 w1 w6 w9 w15]';

%--------------------------------------------------------------------------
%   Points for datum definition
%--------------------------------------------------------------------------          
xy = reshape(X_0(1:8),2,4);

% All points
datum = diag([1 1 1 1]);

%Number of points
p = sum(sum(datum));             

%Centroid
x_c = (1/p)*sum(datum*xy(2,:)');
y_c = (1/p)*sum(datum*xy(1,:)');

%Coordinates reduced to the centroid
y1 = y1-y_c;
y6 = y6-y_c;
y9 = y9-y_c;
y15 = y15-y_c;
x1 = x1-x_c;
x6 = x6-x_c;
x9 = x9-x_c;
x15 = x15-x_c;

%Initial values for unknowns after reduction to the centroid
X_0 = [y1 x1 y6 x6 y9 x9 y15 x15 w1 w6 w9 w15]';

%--------------------------------------------------------------------------     
%Number of observations
no_n=length(L);

%Number of unknowns
no_u=length(X_0);

%Redundancy
r=no_n-no_u+3; %+3 FREE NETWORK...in contrast with exercise 6 where we had control points!!! so no constraint equations!!

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_dist=0.01; %[m] 1cm
%s_dist=0.1; %[m] 10cm APO EKFWNHSH ASKHSHS!!
s_dir=0.001*pi/200;        %Convert to [rad]

s_LL=[s_dist^2*ones(length(dist),1); s_dir^2*ones(length(dir),1)];
S_LL=diag(s_LL);

%Theoretical standard deviation
sigma_0 = 1;

%Cofactor matrix of the observations
Q_LL=1/sigma_0^2*S_LL;

%Weight matrix
P=inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off condition
epsilon=10^-5;
delta=10^-12;
max_x_hat=10^Inf;

%Number of iterations
iteration=0;

%Initialising A
A=zeros(no_n,no_u);

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
    l=L-L_0';

    %Design matrix (14*12)....Prwta oi apostaseis!! to exoume orisei sto vector of observations!!!
    A(1,1)=ds_dy_to(y6,x6,y1,x1);
    A(1,2)=ds_dx_to(y6,x6,y1,x1);
    A(1,3)=ds_dy_from(y6,x6,y1,x1);
    A(1,4)=ds_dx_from(y6,x6,y1,x1);
    
    A(2,1)=ds_dy_to(y9,x9,y1,x1);
    A(2,2)=ds_dx_to(y9,x9,y1,x1);
    A(2,5)=ds_dy_from(y9,x9,y1,x1);
    A(2,6)=ds_dx_from(y9,x9,y1,x1);
    
    A(3,3)=ds_dy_to(y9,x9,y6,x6);
    A(3,4)=ds_dx_to(y9,x9,y6,x6);
    A(3,5)=ds_dy_from(y9,x9,y6,x6);
    A(3,6)=ds_dx_from(y9,x9,y6,x6);
    
    A(4,1)=ds_dy_to(y15,x15,y1,x1);
    A(4,2)=ds_dx_to(y15,x15,y1,x1);
    A(4,7)=ds_dy_from(y15,x15,y1,x1);
    A(4,8)=ds_dx_from(y15,x15,y1,x1);

    A(5,5)=ds_dy_to(y15,x15,y9,x9);
    A(5,6)=ds_dx_to(y15,x15,y9,x9);
    A(5,7)=ds_dy_from(y15,x15,y9,x9);
    A(5,8)=ds_dx_from(y15,x15,y9,x9);
    
    A(6,1)=dr_dy_from(y1,x1,y6,x6);
    A(6,2)=dr_dx_from(y1,x1,y6,x6);
    A(6,3)=dr_dy_to(y1,x1,y6,x6);
    A(6,4)=dr_dx_to(y1,x1,y6,x6);
    A(6,9)=-1;
    
    A(7,1)=dr_dy_from(y1,x1,y15,x15);
    A(7,2)=dr_dx_from(y1,x1,y15,x15);
    A(7,7)=dr_dy_to(y1,x1,y15,x15);
    A(7,8)=dr_dx_to(y1,x1,y15,x15);
    A(7,9)=-1;
    
    A(8,1)=dr_dy_to(y6,x6,y1,x1);
    A(8,2)=dr_dx_to(y6,x6,y1,x1);
    A(8,3)=dr_dy_from(y6,x6,y1,x1);
    A(8,4)=dr_dx_from(y6,x6,y1,x1);
    A(8,10)=-1;
    
    A(9,3)=dr_dy_from(y6,x6,y9,x9);
    A(9,4)=dr_dx_from(y6,x6,y9,x9);
    A(9,5)=dr_dy_to(y6,x6,y9,x9);
    A(9,6)=dr_dx_to(y6,x6,y9,x9);
    A(9,10)=-1;
    
    A(10,5)=dr_dy_from(y9,x9,y15,x15);
    A(10,6)=dr_dx_from(y9,x9,y15,x15);
    A(10,7)=dr_dy_to(y9,x9,y15,x15);
    A(10,8)=dr_dx_to(y9,x9,y15,x15);
    A(10,11)=-1;
    
    A(11,1)=dr_dy_to(y9,x9,y1,x1);
    A(11,2)=dr_dx_to(y9,x9,y1,x1);
    A(11,5)=dr_dy_from(y9,x9,y1,x1);
    A(11,6)=dr_dx_from(y9,x9,y1,x1);
    A(11,11)=-1;
    
    A(12,3)=dr_dy_to(y9,x9,y6,x6);
    A(12,4)=dr_dx_to(y9,x9,y6,x6);
    A(12,5)=dr_dy_from(y9,x9,y6,x6);
    A(12,6)=dr_dx_from(y9,x9,y6,x6);
    A(12,11)=-1;
    
    A(13,1)=dr_dy_to(y15,x15,y1,x1);
    A(13,2)=dr_dx_to(y15,x15,y1,x1);
    A(13,7)=dr_dy_from(y15,x15,y1,x1);
    A(13,8)=dr_dx_from(y15,x15,y1,x1);
    A(13,12)=-1;
    
    A(14,5)=dr_dy_to(y15,x15,y9,x9);
    A(14,6)=dr_dx_to(y15,x15,y9,x9);
    A(14,7)=dr_dy_from(y15,x15,y9,x9);
    A(14,8)=dr_dx_from(y15,x15,y9,x9);
    A(14,12)=-1;
    
    %Normal matrix
     N=A'*P*A; %rank(A)=9!!...The full rank of A it would be 12! So it isn't Maximal!
    %We have 3 columns which are linearly dependent! Rank deficiency=3! So 3 constraints we need!!
    determinant = det(N);
     
    %Matrix B..................... where all the points were contributing to the datum definition
    B=[0 1 0 1 0 1 0 1 0 0 0 0;
       1 0 1 0 1 0 1 0 0 0 0 0;
       -x1 y1 -x6 y6 -x9 y9 -x15 y15 0 0 0 0];
   
    %Extended normal matrix
    N_ext = [N B';B zeros(3,3)];

    %Vector of right hand side of normal equations
    n=A'*P*l;
    
    %Extended vector of right hand side of normal equations
    n_ext = [n;zeros(3,1)];

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx_ext=inv(N_ext);
    
    %Reduced cofactor matrix of the unknowns
    Q_xx = Q_xx_ext(1:no_u,1:no_u);        %Q_11
    
    %Solution of normal equation
    x_hat=Q_xx_ext*n_ext;
    
    %Adjusted unknowns
    X_hat=X_0+x_hat(1:no_u);

    %Update
    X_0=X_hat;

    y1=X_0(1);
    x1=X_0(2);
    y6=X_0(3);
    x6=X_0(4);
    y9=X_0(5);
    x9=X_0(6);
    y15=X_0(7);
    x15=X_0(8);
    w1=X_0(9);
    w6=X_0(10);
    w9=X_0(11);
    w15=X_0(12);

    %Check 1
    max_x_hat=max(abs(x_hat));
	
	%Vector of residuals
    v = A*x_hat(1:no_u)-l;
    v_gon = v(6:14,1)*200/pi;         %Convert to [gon]

    %Objective function
    vTPv = v'*P*v;

    %Vector of adjusted observations
    L_hat = L+v;
    L_hat_gon = L_hat(6:14)*200/pi;   %Convert to [gon]

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
    iteration=iteration+1;

end

% Final Check
if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

%Empirical reference standard deviation
s_0=sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat=s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));
s_X_gon = s_X(9:12,1)*200/pi;     %Convert to [gon]

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));
s_L_hat_gon = s_L_hat(6:14)*200/pi;    %Convert to [gon]

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));
s_v_gon = s_v(6:14,1)*200/pi;         %Convert to [gon]

%De-reduced unknowns
y1 = X_0(1)+y_c;
x1 = X_0(2)+x_c;
y6 = X_0(3)+y_c;
x6 = X_0(4)+x_c;
y9 = X_0(5)+y_c;
x9 = X_0(6)+x_c;
y15 = X_0(7)+y_c;
x15 = X_0(8)+x_c;

X_final = [y1 x1 y6 x6 y9 x9 y15 x15 w1 w6 w9 w15]';

%Convert to [gon] and check the quadrants
gon = X_final(9:12,1)*200/pi;
gon(1) = gon(1)+400;
gon(2) = gon(2)+400;
gon(4) = gon(4)+400;




%% Global test, s=95%...a=5%
Tx2 = (r*s_0^2)/sigma_0^2;

%Two-sided problem
%1-a/2 = 1-0.05/2=0.975
%a/2 = 0.05/2 = 0.025

tx2u = chi2inv(0.975,r); %INVERSE FUNCTION!!
tx2l = chi2inv(0.025,r);

if tx2l<Tx2 && Tx2<tx2u
  disp('The test fails to reject the Ho.')
else
  disp('The test rejects the Ho.')
end  

%--------------------------------------------------------------------------
% Internal and external reliability parameters
%--------------------------------------------------------------------------
% Parameters for internal
% Redundancy numbers
EV = 100*(diag(Q_vv*P)); %Eveyrthing is fine since all redundancy numbers are between 10-70%!!!

% Standardised residuals
sigma_v = sigma_0*sqrt(diag(Q_vv));
NV = abs(v) ./ sigma_v; %No blunder distinguishable since NV<2.5

% Potential magnitude of a blunder
GF = -v./(diag(Q_vv*P)); %The last values are very small because are the directions in rad
GF_gon = GF(6:14,1)*200/pi; %convert to [gon]

% Lower boundary value for blunders
GRZW = ones(no_n,1);
for i=1:no_n
  GRZW(i,1) = sigma_0*4.13/(sqrt(EV(i,1)*P(i,i)));
end
GRZW_gon = GRZW(6:14,1)*200/pi; %convert to [gon]

%Parameters for external reliability
P_diag = diag(P);
r_w = ones(9,1); %since 9 are directions observations

%1st data set 
%ston PARONOMASTI sum of the weights of all directions of one data set
r_w(1) = P_diag(6,1)/(P_diag(6,1)+P_diag(7,1)); %sum of weights of one data set FROM 1 TO 6, FROM 1 TO 15!
r_w(2) = P_diag(7,1)/(P_diag(6,1)+P_diag(7,1)); 

%2nd data set
r_w(3) = P_diag(8,1)/(P_diag(8,1)+P_diag(9,1));
r_w(4) = P_diag(9,1)/(P_diag(8,1)+P_diag(9,1));
%3rd data set
r_w(5) = P_diag(10,1)/(P_diag(10,1)+P_diag(11,1)+P_diag(12,1));
r_w(6) = P_diag(11,1)/(P_diag(10,1)+P_diag(11,1)+P_diag(12,1));
r_w(7) = P_diag(12,1)/(P_diag(10,1)+P_diag(11,1)+P_diag(12,1));
%4th data set
r_w(8) = P_diag(13,1)/(P_diag(13,1)+P_diag(14,1));
r_w(9) = P_diag(14,1)/(P_diag(13,1)+P_diag(14,1));

%Distances because we have to do s/p  %organise it according to table1 dist
dd(1) = dist(1,3); %the first distance observation...3 since is in the 3rd column
dd(2) = dist(4,3);
dd(3) = dist(1,3); %FROM 1 TO 6 is the same with the DISTANCE FROM 6 TO 1
dd(4) = dist(3,3);
dd(5) = dist(5,3);
dd(6) = dist(2,3);
dd(7) = dist(3,3);
dd(8) = dist(4,3);
dd(9) = dist(5,3);

%Impact of the boundary value on the coordinates of the corresponding
%points
EGK = ones(no_n,1);

%For distances!!!
for i=1:5
  EGK(i,1)=(1-EV(i,1))*GRZW(i,1);
end

%for directions!!!
for i=6:14
  EGK(i,1)=(1-EV(i,1)-r_w(i-5,1))*GRZW(i,1)*dd(i-5); %i-5 because the iteration begins from 6:14
end

%Impact of a potential blunder on a point corresponding to the measurement

EP = ones(14,1);

%For distances!!!
for i=1:5
  EP(i,1)=(1-EV(i,1))*GF(i,1);
end

%for directions!!!
%i-5 because the iteration begins from 6:14
for i=6:14
  EP(i,1)=(1-EV(i,1)-r_w(i-5,1))*GF(i,1)*dd(i-5);
end

%% Condition for observation removal.

% Initialize a flag variable
displayed = false;

for i = 1:14
    if EV(i) >= 30
        if NV(i) < 4 && abs(EP(i)) < 0.002
            if ~displayed
                disp('The result satisfies the requirements for precision and reliability');
                displayed = true; % Set the flag to true
            end
        elseif NV(i) > 4 && abs(EP(i)) < 0.002
            if ~displayed
                disp('The observation does not have to be removed since it will only have a minor impact');
                displayed = true; % Set the flag to true
            end
        elseif NV(i) > 4 && abs(EP(i)) > 0.002
            max_NV_index = find(NV == max(NV));
            fprintf('The %d-th observation needs to be removed and re-adjusted\n', max_NV_index);
        end
    elseif EV(i) < 10 && NV(i) < 4 && abs(EP(i)) < 0.002
        fprintf('We have to introduce more observations; although the observation at index %d meets the criteria, EV is less than 10%%\n', i);
    end
end








