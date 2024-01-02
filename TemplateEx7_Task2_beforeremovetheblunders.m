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
dist = load('Distances_Task2.txt'); %different network than task 1
dir = load('Directions_Task2.txt');

%Vector of observations
L = [dist(:,3); dir(:,3)*pi/200];    %Convert to [rad]

%Gauss-Krueger coordinates for new points [m]
y1000 = 4590337.390; %Easting [m]
x1000 = 5820823.642;
y2000 = 4589967.526;
x2000 = 5820806.067;
y3000 = 4590078.021;
x3000 = 5820681.807;

y100 = 4590159.8;
x100 = 5820727.3;
y101 = 4589800.1;
x101 = 5820857.9;
y102 = 4590163.2;
x102 = 5820848.7;
y103 = 4589956.9;
x103 = 5820700.4;

%Initial values for orientation unknowns
w1000 = 0;
w2000 = 0;
w3000 = 0;
w100 = 0;
w101 = 0;
w102 = 0;
w103 = 0;

%Initial values for unknowns
X_0 = [y1000 x1000 y2000 x2000 y3000 x3000 y100 x100 y101 x101 y102 x102 y103 x103 w1000 w2000 w3000 w100 w101 w102 w103]';

%--------------------------------------------------------------------------
%   Points for datum definition
%--------------------------------------------------------------------------          
xy = reshape(X_0(1:14),2,7); %1st row includes y-Easting

% All points
datum = diag([1 1 1 0 0 0 0]); %points of table 3 [1000, 2000, 3000] are contributing to the datum definition

%Number of points
p = sum(sum(datum));             

%Centroid
x_c = (1/p)*sum(datum*xy(2,:)');
y_c = (1/p)*sum(datum*xy(1,:)');

%Coordinates reduced to the centroid
y1000 = y1000-y_c;
y2000 = y2000-y_c;
y3000 = y3000-y_c;
y100 = y100-y_c;
y101 = y101-y_c;
y102 = y102-y_c;
y103 = y103-y_c;
x1000 = x1000-x_c;
x2000 = x2000-x_c;
x3000 = x3000-x_c;
x100 = x100-x_c;
x101 = x101-x_c;
x102 = x102-x_c;
x103 = x103-x_c;

%Initial values for unknowns after reduction to the centroid
X_0 = [y1000 x1000 y2000 x2000 y3000 x3000 y100 x100 y101 x101 y102 x102 y103 x103 w1000 w2000 w3000 w100 w101 w102 w103]';

%--------------------------------------------------------------------------     
%Number of observations
no_n=length(L);

%Number of unknowns
no_u=length(X_0);

%Redundancy
r=no_n-no_u+3; %18

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_dist=0.001;   % [m] if we need to fail to reject we have to use 1cm!!
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
    L_0(1)=dis(y1000,x1000,y100,x100);
    L_0(2)=dis(y1000,x1000,y102,x102);
    L_0(3)=dis(y2000,x2000,y103,x103); %this is a blunder!!!
    L_0(4)=dis(y2000,x2000,y101,x101);
    L_0(5)=dis(y3000,x3000,y100,x100);
    L_0(6)=dis(y3000,x3000,y103,x103);
    L_0(7)=dis(y100,x100,y1000,x1000);
    L_0(8)=dis(y100,x100,y102,x102);
    L_0(9)=dis(y100,x100,y2000,x2000);
    L_0(10)=dis(y100,x100,y3000,x3000);
    L_0(11)=dis(y101,x101,y103,x103);
    L_0(12)=dis(y101,x101,y2000,x2000);
    L_0(13)=dis(y102,x102,y1000,x1000);
    L_0(14)=dis(y102,x102,y100,x100);
    L_0(15)=dis(y102,x102,y2000,x2000);
    L_0(16)=dis(y103,x103,y3000,x3000);
    L_0(17)=dis(y103,x103,y2000,x2000);
    L_0(18)=dis(y103,x103,y101,x101);
    
    %Vector of directions
    L_0(19)=direction(y1000,x1000,y100,x100,w1000);
    L_0(20)=direction(y1000,x1000,y102,x102,w1000);
    L_0(21)=direction(y2000,x2000,y103,x103,w2000);
    L_0(22)=direction(y2000,x2000,y101,x101,w2000);
    L_0(23)=direction(y3000,x3000,y100,x100,w3000);
    L_0(24)=direction(y3000,x3000,y103,x103,w3000);
    L_0(25)=direction(y100,x100,y1000,x1000,w100);
    L_0(26)=direction(y100,x100,y102,x102,w100); %this is also a blunder!!!
    L_0(27)=direction(y100,x100,y2000,x2000,w100);
    L_0(28)=direction(y100,x100,y3000,x3000,w100);
    L_0(29)=direction(y101,x101,y103,x103,w101);
    L_0(30)=direction(y101,x101,y2000,x2000,w101);
    L_0(31)=direction(y102,x102,y1000,x1000,w102);
    L_0(32)=direction(y102,x102,y100,x100,w102);
    L_0(33)=direction(y102,x102,y2000,x2000,w102);
    L_0(34)=direction(y103,x103,y3000,x3000,w103);
    L_0(35)=direction(y103,x103,y2000,x2000,w103);
    L_0(36)=direction(y103,x103,y101,x101,w103);
    
    %Vector of reduced observations
    l=L-L_0'; %(36*1)

    %Design matrix
    A(1,1) = ds_dy_from(y1000,x1000,y100,x100); %From 1000 TO 100
    A(1,2) = ds_dx_from(y1000,x1000,y100,x100);
    A(1,7) = ds_dy_to(y1000,x1000,y100,x100);
    A(1,8) = ds_dx_to(y1000,x1000,y100,x100);
    
    A(2,1) = ds_dy_from(y1000,x1000,y102,x102);
    A(2,2) = ds_dx_from(y1000,x1000,y102,x102);
    A(2,11) = ds_dy_to(y1000,x1000,y102,x102);
    A(2,12) = ds_dx_to(y1000,x1000,y102,x102);
    
    A(3,3) = ds_dy_from(y2000,x2000,y103,x103);
    A(3,4) = ds_dx_from(y2000,x2000,y103,x103);
    A(3,13) = ds_dy_to(y2000,x2000,y103,x103);
    A(3,14) = ds_dx_to(y2000,x2000,y103,x103);
    
    A(4,3) = ds_dy_from(y2000,x2000,y101,x101);
    A(4,4) = ds_dx_from(y2000,x2000,y101,x101);
    A(4,9) = ds_dy_to(y2000,x2000,y101,x101);
    A(4,10) = ds_dx_to(y2000,x2000,y101,x101);
    
    A(5,5) = ds_dy_from(y3000,x3000,y100,x100);
    A(5,6) = ds_dx_from(y3000,x3000,y100,x100);
    A(5,7) = ds_dy_to(y3000,x3000,y100,x100);
    A(5,8) = ds_dx_to(y3000,x3000,y100,x100);
    
    A(6,5) = ds_dy_from(y3000,x3000,y103,x103);
    A(6,6) = ds_dx_from(y3000,x3000,y103,x103);
    A(6,13) = ds_dy_to(y3000,x3000,y103,x103);
    A(6,14) = ds_dx_to(y3000,x3000,y103,x103);
    
    A(7,1) = ds_dy_to(y100,x100,y1000,x1000);
    A(7,2) = ds_dx_to(y100,x100,y1000,x1000);
    A(7,7) = ds_dy_from(y100,x100,y1000,x1000);
    A(7,8) = ds_dx_from(y100,x100,y1000,x1000);
    
    A(8,7) = ds_dy_from(y100,x100,y102,x102);
    A(8,8) = ds_dx_from(y100,x100,y102,x102);
    A(8,11) = ds_dy_to(y100,x100,y102,x102);
    A(8,12) = ds_dx_to(y100,x100,y102,x102);
    
    A(9,3) = ds_dy_to(y100,x100,y2000,x2000);
    A(9,4) = ds_dx_to(y100,x100,y2000,x2000);
    A(9,7) = ds_dy_from(y100,x100,y2000,x2000);
    A(9,8) = ds_dx_from(y100,x100,y2000,x2000);
    
    A(10,5) = ds_dy_to(y100,x100,y3000,x3000); %10th observation of DISTANCES....FROM 100 TO 3000
    A(10,6) = ds_dx_to(y100,x100,y3000,x3000);
    A(10,7) = ds_dy_from(y100,x100,y3000,x3000);
    A(10,8) = ds_dx_from(y100,x100,y3000,x3000);
    
    A(11,9) = ds_dy_from(y101,x101,y103,x103);
    A(11,10) = ds_dx_from(y101,x101,y103,x103);
    A(11,13) = ds_dy_to(y101,x101,y103,x103);
    A(11,14) = ds_dx_to(y101,x101,y103,x103);
    
    A(12,3) = ds_dy_to(y101,x101,y2000,x2000);
    A(12,4) = ds_dx_to(y101,x101,y2000,x2000);
    A(12,9) = ds_dy_from(y101,x101,y2000,x2000);
    A(12,10) = ds_dx_from(y101,x101,y2000,x2000);
    
    A(13,1) = ds_dy_to(y102,x102,y1000,x1000);
    A(13,2) = ds_dx_to(y102,x102,y1000,x1000);
    A(13,11) = ds_dy_from(y102,x102,y1000,x1000);
    A(13,12) = ds_dx_from(y102,x102,y1000,x1000);
    
    A(14,7) = ds_dy_to(y102,x102,y100,x100);
    A(14,8) = ds_dx_to(y102,x102,y100,x100);
    A(14,11) = ds_dy_from(y102,x102,y100,x100);
    A(14,12) = ds_dx_from(y102,x102,y100,x100);
    
    A(15,3) = ds_dy_to(y102,x102,y2000,x2000);
    A(15,4) = ds_dx_to(y102,x102,y2000,x2000);
    A(15,11) = ds_dy_from(y102,x102,y2000,x2000);
    A(15,12) = ds_dx_from(y102,x102,y2000,x2000);
    
    A(16,5) = ds_dy_to(y103,x103,y3000,x3000); %DISTANCE FROM 103 TO 3000
    A(16,6) = ds_dx_to(y103,x103,y3000,x3000);
    A(16,13) = ds_dy_from(y103,x103,y3000,x3000);
    A(16,14) = ds_dx_from(y103,x103,y3000,x3000);
    
    A(17,3) = ds_dy_to(y103,x103,y2000,x2000);
    A(17,4) = ds_dx_to(y103,x103,y2000,x2000);
    A(17,13) = ds_dy_from(y103,x103,y2000,x2000);
    A(17,14) = ds_dx_from(y103,x103,y2000,x2000);
    
    A(18,9) = ds_dy_to(y103,x103,y101,x101); %lAST Distance observation FROM 103 TO 101 
    A(18,10) = ds_dx_to(y103,x103,y101,x101);
    A(18,13) = ds_dy_from(y103,x103,y101,x101);
    A(18,14) = ds_dx_from(y103,x103,y101,x101);
    
    A(19,1) = dr_dy_from(y1000,x1000,y100,x100);
    A(19,2) = dr_dx_from(y1000,x1000,y100,x100);
    A(19,7) = dr_dy_to(y1000,x1000,y100,x100);
    A(19,8) = dr_dx_to(y1000,x1000,y100,x100);
    A(19,15) = -1;
    
    A(20,1) = dr_dy_from(y1000,x1000,y102,x102);
    A(20,2) = dr_dx_from(y1000,x1000,y102,x102);
    A(20,11) = dr_dy_to(y1000,x1000,y102,x102);
    A(20,12) = dr_dx_to(y1000,x1000,y102,x102);
    A(20,15) = -1;
    
    A(21,3) = dr_dy_from(y2000,x2000,y103,x103);
    A(21,4) = dr_dx_from(y2000,x2000,y103,x103);
    A(21,13) = dr_dy_to(y2000,x2000,y103,x103);
    A(21,14) = dr_dx_to(y2000,x2000,y103,x103);
    A(21,16) = -1;
    
    A(22,3) = dr_dy_from(y2000,x2000,y101,x101);
    A(22,4) = dr_dx_from(y2000,x2000,y101,x101);
    A(22,9) = dr_dy_to(y2000,x2000,y101,x101);
    A(22,10) = dr_dx_to(y2000,x2000,y101,x101);
    A(22,16) = -1;
    
    A(23,5) = dr_dy_from(y3000,x3000,y100,x100);
    A(23,6) = dr_dx_from(y3000,x3000,y100,x100);
    A(23,7) = dr_dy_to(y3000,x3000,y100,x100);
    A(23,8) = dr_dx_to(y3000,x3000,y100,x100);
    A(23,17) = -1;
    
    A(24,5) = dr_dy_from(y3000,x3000,y103,x103);
    A(24,6) = dr_dx_from(y3000,x3000,y103,x103);
    A(24,13) = dr_dy_to(y3000,x3000,y103,x103);
    A(24,14) = dr_dx_to(y3000,x3000,y103,x103);
    A(24,17) = -1;
    
    A(25,1) = dr_dy_to(y100,x100,y1000,x1000);
    A(25,2) = dr_dx_to(y100,x100,y1000,x1000);
    A(25,7) = dr_dy_from(y100,x100,y1000,x1000);
    A(25,8) = dr_dx_from(y100,x100,y1000,x1000);
    A(25,18) = -1;
    
    A(26,7) = dr_dy_from(y100,x100,y102,x102);
    A(26,8) = dr_dx_from(y100,x100,y102,x102);
    A(26,11) = dr_dy_to(y100,x100,y102,x102);
    A(26,12) = dr_dx_to(y100,x100,y102,x102);
    A(26,18) = -1;
    
    A(27,3) = dr_dy_to(y100,x100,y2000,x2000);
    A(27,4) = dr_dx_to(y100,x100,y2000,x2000);
    A(27,7) = dr_dy_from(y100,x100,y2000,x2000);
    A(27,8) = dr_dx_from(y100,x100,y2000,x2000);
    A(27,18) = -1;
    
    A(28,5) = dr_dy_to(y100,x100,y3000,x3000);
    A(28,6) = dr_dx_to(y100,x100,y3000,x3000);
    A(28,7) = dr_dy_from(y100,x100,y3000,x3000);
    A(28,8) = dr_dx_from(y100,x100,y3000,x3000);
    A(28,18) = -1;
    
    A(29,9) = dr_dy_from(y101,x101,y103,x103);
    A(29,10) = dr_dx_from(y101,x101,y103,x103);
    A(29,13) = dr_dy_to(y101,x101,y103,x103);
    A(29,14) = dr_dx_to(y101,x101,y103,x103);
    A(29,19) = -1;
    
    A(30,3) = dr_dy_to(y101,x101,y2000,x2000);
    A(30,4) = dr_dx_to(y101,x101,y2000,x2000);
    A(30,9) = dr_dy_from(y101,x101,y2000,x2000);
    A(30,10) = dr_dx_from(y101,x101,y2000,x2000);
    A(30,19) = -1;
    
    A(31,1) = dr_dy_to(y102,x102,y1000,x1000);
    A(31,2) = dr_dx_to(y102,x102,y1000,x1000);
    A(31,11) = dr_dy_from(y102,x102,y1000,x1000);
    A(31,12) = dr_dx_from(y102,x102,y1000,x1000);
    A(31,20) = -1;
    
    A(32,7) = dr_dy_to(y102,x102,y100,x100);
    A(32,8) = dr_dx_to(y102,x102,y100,x100);
    A(32,11) = dr_dy_from(y102,x102,y100,x100);
    A(32,12) = dr_dx_from(y102,x102,y100,x100);
    A(32,20) = -1;
    
    A(33,3) = dr_dy_to(y102,x102,y2000,x2000);
    A(33,4) = dr_dx_to(y102,x102,y2000,x2000);
    A(33,11) = dr_dy_from(y102,x102,y2000,x2000);
    A(33,12) = dr_dx_from(y102,x102,y2000,x2000);
    A(33,20) = -1;
    
    A(34,5) = dr_dy_to(y103,x103,y3000,x3000);
    A(34,6) = dr_dx_to(y103,x103,y3000,x3000);
    A(34,13) = dr_dy_from(y103,x103,y3000,x3000);
    A(34,14) = dr_dx_from(y103,x103,y3000,x3000);
    A(34,21) = -1;
    
    A(35,3) = dr_dy_to(y103,x103,y2000,x2000);
    A(35,4) = dr_dx_to(y103,x103,y2000,x2000);
    A(35,13) = dr_dy_from(y103,x103,y2000,x2000);
    A(35,14) = dr_dx_from(y103,x103,y2000,x2000);
    A(35,21) = -1;
    
    A(36,9) = dr_dy_to(y103,x103,y101,x101); %LAST OBSERVATION FROM 103 TO 101
    A(36,10) = dr_dx_to(y103,x103,y101,x101);
    A(36,13) = dr_dy_from(y103,x103,y101,x101);
    A(36,14) = dr_dx_from(y103,x103,y101,x101);
    A(36,21) = -1; %orientation parameter EVER in the point FROM!!
    
    
    %Normal matrix
     N=A'*P*A;
     rankA = rank(A); %rank(A)=18!!...The full rank of A it would be 21! So it isn't Maximal!
    %We have 3 columns which are linearly dependent! Rank deficiency=3! So 3 constraints we need!!
     
    %Matrix B ....all points in Table 3 [1000,2000,3000] are contributing to the datum definition
    B = [0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         -x1000 y1000 -x2000 y2000 -x3000 y3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

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

    y1000 = X_0(1);
    x1000 = X_0(2);
    y2000 = X_0(3);
    x2000 = X_0(4);
    y3000 = X_0(5);
    x3000 = X_0(6);
    y100 = X_0(7);
    x100 = X_0(8);
    y101 = X_0(9);
    x101 = X_0(10);
    y102 = X_0(11);
    x102 = X_0(12);
    y103 = X_0(13);
    x103 = X_0(14);
    w1000 = X_0(15);
    w2000 = X_0(16);
    w3000 = X_0(17);
    w100 = X_0(18);
    w101 = X_0(19);
    w102 = X_0(20);
    w103 = X_0(21);

    %Check 1
    max_x_hat=max(abs(x_hat));
	
	%Vector of residuals
    v = A*x_hat(1:no_u)-l;
    v_gon = v(19:36,1)*200/pi;         %Convert to [gon]

    %Objective function
    vTPv = v'*P*v;

    %Vector of adjusted observations
    L_hat = L+v;
    L_hat_gon = L_hat(19:36)*200/pi;   %Convert to [gon]

    % distances
    Phi_X_hat(1) = dis(y1000,x1000,y100,x100);
    Phi_X_hat(2)=dis(y1000,x1000,y102,x102);
    Phi_X_hat(3)=dis(y2000,x2000,y103,x103); %it is blunder
    Phi_X_hat(4)=dis(y2000,x2000,y101,x101);
    Phi_X_hat(5)=dis(y3000,x3000,y100,x100);
    Phi_X_hat(6)=dis(y3000,x3000,y103,x103);
    Phi_X_hat(7)=dis(y100,x100,y1000,x1000);
    Phi_X_hat(8)=dis(y100,x100,y102,x102);
    Phi_X_hat(9)=dis(y100,x100,y2000,x2000);
    Phi_X_hat(10)=dis(y100,x100,y3000,x3000);
    Phi_X_hat(11)=dis(y101,x101,y103,x103);
    Phi_X_hat(12)=dis(y101,x101,y2000,x2000);
    Phi_X_hat(13)=dis(y102,x102,y1000,x1000);
    Phi_X_hat(14)=dis(y102,x102,y100,x100);
    Phi_X_hat(15)=dis(y102,x102,y2000,x2000);
    Phi_X_hat(16)=dis(y103,x103,y3000,x3000);
    Phi_X_hat(17)=dis(y103,x103,y2000,x2000);
    Phi_X_hat(18)=dis(y103,x103,y101,x101);
	
    % directions
    Phi_X_hat(19) = direction(y1000,x1000,y100,x100,w1000);
    Phi_X_hat(20)=direction(y1000,x1000,y102,x102,w1000);
    Phi_X_hat(21)=direction(y2000,x2000,y103,x103,w2000);
    Phi_X_hat(22)=direction(y2000,x2000,y101,x101,w2000);
    Phi_X_hat(23)=direction(y3000,x3000,y100,x100,w3000);
    Phi_X_hat(24)=direction(y3000,x3000,y103,x103,w3000);
    Phi_X_hat(25)=direction(y100,x100,y1000,x1000,w100);
    Phi_X_hat(26)=direction(y100,x100,y102,x102,w100); %it is also a blunder!!! it will be removed!!
    Phi_X_hat(27)=direction(y100,x100,y2000,x2000,w100);
    Phi_X_hat(28)=direction(y100,x100,y3000,x3000,w100);
    Phi_X_hat(29)=direction(y101,x101,y103,x103,w101);
    Phi_X_hat(30)=direction(y101,x101,y2000,x2000,w101);
    Phi_X_hat(31)=direction(y102,x102,y1000,x1000,w102);
    Phi_X_hat(32)=direction(y102,x102,y100,x100,w102);
    Phi_X_hat(33)=direction(y102,x102,y2000,x2000,w102);
    Phi_X_hat(34)=direction(y103,x103,y3000,x3000,w103);
    Phi_X_hat(35)=direction(y103,x103,y2000,x2000,w103);
    Phi_X_hat(36)=direction(y103,x103,y101,x101,w103);
  
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

%Empirical reference standard deviation.....
s_0=sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat=s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));
s_X_gon = s_X(15:21,1)*200/pi;     %Convert to [gon]

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));
s_L_hat_gon = s_L_hat(19:36)*200/pi;    %Convert to [gon]

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));
s_v_gon = s_v(19:36,1)*200/pi;         %Convert to [gon]

%De-reduced unknowns
y1000 = X_0(1)+y_c;
x1000 = X_0(2)+x_c;
y2000 = X_0(3)+y_c;
x2000 = X_0(4)+x_c;
y3000 = X_0(5)+y_c;
x3000 = X_0(6)+x_c;
y100 = X_0(7)+y_c;
x100 = X_0(8)+x_c;
y101 = X_0(9)+y_c;
x101 = X_0(10)+x_c;
y102 = X_0(11)+y_c;
x102 = X_0(12)+x_c;
y103 = X_0(13)+y_c;
x103 = X_0(14)+x_c;

X_final =  [y1000 x1000 y2000 x2000 y3000 x3000 y100 x100 y101 x101 y102 x102 y103 x103 w1000 w2000 w3000 w100 w101 w102 w103]';

%Convert to [gon] and check the quadrants
gon = X_final(15:21,1)*200/pi;
gon(1) = gon(1)+400;
gon(2) = gon(2)+400;
gon(3) = gon(3)+400;
gon(4) = gon(4)+400;
gon(5) = gon(5)+400;
gon(6) = gon(6)+400;
gon(7) = gon(7)+400;




%% Global test, S = 95%, a = 5%.............................................................................
Tx2 = (r*s_0^2)/sigma_0^2;

%Two-sided problem 
%1-a/2 = 1-0.05/2 = 0.975
%a/2 = 0.05/2 = 0.025

tx2u = chi2inv(0.975, r); %INVERSE FUNCTION!! X2
tx2l = chi2inv(0.025, r);

if tx2l<Tx2 && Tx2<tx2u
  disp('The test fails to reject the H0.')
else
  disp('The test rejects the H0.') %Initially the test rejects the null hypothesis!! So Significantly different!!!
end

%--------------------------------------------------------------------------
% Internal and external reliability parameters
%--------------------------------------------------------------------------
% Parameters for internal
% Redundancy numbers
EV = 100*(diag(Q_vv*P)); %Eveyrthing is fine if all redundancy numbers are between 30-70%!!!

% Standardised residuals
sigma_v = sigma_0*sqrt(diag(Q_vv));
NV = abs(v)./sigma_v;

% Potential magnitude of a blunder
GF = -v./(diag(Q_vv*P));
GF_gon = GF(19:36,1)*200/pi; %convert to [gon]

% Lower boundary value for blunders
GRZW = ones(no_n,1);
for i=1:no_n
  GRZW(i,1) = sigma_0*4.13/(sqrt(EV(i,1)*P(i,i)));  %4.13-constant
end  % EV(i,1) is the same as ri x 100%

GRZW_gon = GRZW(19:36,1)*200/pi; %convert to [gon]

%Parameters for external reliability
P_diag = diag(P);
r_w = ones(18,1);  %18 directions

%ston PARONOMASTI sum of the weights of all directions of one data set
%dataset 1 - from point 1000  % In L, order 19-36 are directions
r_w(1) = P_diag(19, 1)/(P_diag(19, 1)+P_diag(20, 1)); %sum of weights of one data set FROM 1000 TO 100, FROM 1000 TO 102!!!
r_w(2) = P_diag(20, 1)/(P_diag(19, 1)+P_diag(20, 1));

%dataset 2 - from point 2000
r_w(3) = P_diag(21, 1)/(P_diag(21, 1)+P_diag(22, 1));
r_w(4) = P_diag(22, 1)/(P_diag(21, 1)+P_diag(22, 1));

%dataset 3 - from point 3000
r_w(5) = P_diag(23, 1)/(P_diag(23, 1)+P_diag(24, 1));
r_w(6) = P_diag(24, 1)/(P_diag(23, 1)+P_diag(24, 1));

%dataset 4 - from point 100
r_w(7) = P_diag(25, 1)/(P_diag(25, 1)+P_diag(26, 1)+P_diag(27, 1)+P_diag(28, 1));
r_w(8) = P_diag(26, 1)/(P_diag(25, 1)+P_diag(26, 1)+P_diag(27, 1)+P_diag(28, 1)); %It's deleted later since 100-102 is blunder!!!
r_w(9) = P_diag(27, 1)/(P_diag(25, 1)+P_diag(26, 1)+P_diag(27, 1)+P_diag(28, 1));
r_w(10) = P_diag(28, 1)/(P_diag(25, 1)+P_diag(26, 1)+P_diag(27, 1)+P_diag(28, 1));

%dataset 5 - from point 101
r_w(11) = P_diag(29, 1)/(P_diag(29, 1)+P_diag(30, 1));
r_w(12) = P_diag(30, 1)/(P_diag(29, 1)+P_diag(30, 1));

%dataset 6 - from point 102
r_w(13) = P_diag(31, 1)/(P_diag(31, 1)+P_diag(32, 1)+P_diag(33, 1));
r_w(14) = P_diag(32, 1)/(P_diag(31, 1)+P_diag(32, 1)+P_diag(33, 1));
r_w(15) = P_diag(33, 1)/(P_diag(31, 1)+P_diag(32, 1)+P_diag(33, 1));

%dataset 7 - from point 103
r_w(16) = P_diag(34, 1)/(P_diag(34, 1)+P_diag(35, 1)+P_diag(36, 1));
r_w(17) = P_diag(35, 1)/(P_diag(34, 1)+P_diag(35, 1)+P_diag(36, 1));
r_w(18) = P_diag(36, 1)/(P_diag(34, 1)+P_diag(35, 1)+P_diag(36, 1));

%Distances
%Distances
for i = 1:18
  dd(i) = dist(i, 3);
end

%Impact of the boundary value on the coordinates of the corresponding
%points
EGK = ones(36,1);
%for distances
for i=1:18
  EGK(i,1) = (1-EV(i,1))*GRZW(i,1);
end

%for directions
for i=19:36
  EGK(i,1) = (1-EV(i,1) - r_w(i-18,1))*GRZW(i,1)*dd(i-18);
end  

%Impact of a potential blunder on a point corresponding to the measurement
EP = ones(36,1);

%for distances
for i=1:18
  EP(i,1) = (1-EV(i,1))*GF(i,1);
end

%for directions
for i=19:36
  EP(i,1) = (1 - EV(i,1) - r_w(i-18,1))*GF(i,1)*dd(i-18);
end  









