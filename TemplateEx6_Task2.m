%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION
%    Network adjustment with stochastic datum     
% 
%   Author         : Anastasia Pasioti
%   Version        : June 7, 2017
%   Last changes   : June 4, 2018
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format long g;

%Load all files
dir = load('Directions_Task2.txt');
dist = load('Distances_Task2.txt');
fixedpoint = load('ControlPoints_Task2.txt'); %Gauss-Krueger coordinates of the control points...table3
newpoint = load('NewPoints_Task2.txt'); %Approximate values for the coordinates of the new points...table4

%Coordinates of the control points
y1000 = fixedpoint(1,2);
x1000 = fixedpoint(1,3);
y2000 = fixedpoint(2,2);
x2000 = fixedpoint(2,3);
y3000 = fixedpoint(3,2);
x3000 = fixedpoint(3,3);
%Initial values for the new points
y100 = newpoint(1,2);
x100 = newpoint(1,3);
y101 = newpoint(2,2);
x101 = newpoint(2,3);
y102 = newpoint(3,2);
x102 = newpoint(3,3);
y103 = newpoint(4,2);
x103 = newpoint(4,3);
%Initial values for the orientation unknowns
w1000 = 0; %Linear terms so initial value=0
w2000 = 0;
w3000 = 0;
w100 = 0;
w101 = 0;
w102 = 0;
w103 = 0;

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
%Vector of Observations
L = [dist(:,3); dir(:,3)*pi/200; y1000; x1000; y2000; x2000; y3000; x3000]; %I will include also the coordinates of the control points!!
%Control points will be introcuded as unknowns and as observations

%Initial values for unknowns
X_0 = [y1000 x1000 y2000 x2000 y3000 x3000 y100 x100 y101 x101 y102 x102 y103 x103 w1000 w2000 w3000 w100 w101 w102 w103]';

%Number of observations
no_n = length(L);

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n - no_u;   %WE DO NOT HAVE CONSTRAINT EQUATIONS!!!

xy = [y1000 x1000 y2000 x2000 y3000 x3000]; %We write this because of (length(xy),1)...otherwise (6,1)
%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_dist = 0.001;   %[m]
s_dir = 0.001*pi/200;   %[gon] -> [rad]
s_xy = 0.01;    %[m]Std for the coordinates of control points!! Assumption for our adjustment

%Variance-Covariance Matrix of the observations
S_LL = diag([s_dist^2*ones(length(dist),1); s_dir^2*ones(length(dir),1); s_xy^2*ones(length(xy),1)]);

%Theoretical standard deviation
sigma_0 = 1;

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off condition
epsilon = 10^-5;
delta = 10^-9;
max_x_hat = 10^Inf; %Of course very large number!!!

%Number of iterations
iteration = 0;

%Initialization A
A = zeros(no_n,no_u);

%Iteration
while max_x_hat>epsilon || Check2>delta
%................Functional Model..............
   %Vector of distances.....18 measurments of distances!!
    L_0(1)=dis(y1000,x1000,y100,x100); %From 1000 to 100
    L_0(2)=dis(y1000,x1000,y102,x102); %From 1000 to 102
    L_0(3)=dis(y2000,x2000,y103,x103); %From 2000 to 103
    L_0(4)=dis(y2000,x2000,y101,x101); %From 2000 to 101
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
    
    %Vector of directions.....18 measurments of directions!!
    L_0(19)=direction(y1000,x1000,y100,x100,w1000); %From 1000 to 100
    L_0(20)=direction(y1000,x1000,y102,x102,w1000);
    L_0(21)=direction(y2000,x2000,y103,x103,w2000);
    L_0(22)=direction(y2000,x2000,y101,x101,w2000);
    L_0(23)=direction(y3000,x3000,y100,x100,w3000);
    L_0(24)=direction(y3000,x3000,y103,x103,w3000);
    L_0(25)=direction(y100,x100,y1000,x1000,w100);
    L_0(26)=direction(y100,x100,y102,x102,w100);
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
    
    %Observed unknowns....coordinates of control points
    L_0(37) = y1000;
    L_0(38) = x1000;
    L_0(39) = y2000;
    L_0(40) = x2000;
    L_0(41) = y3000;
    L_0(42) = x3000;
    
    %Vector of reduced observations
    l = L-L_0';

    %Design matrix
    A(1,1) = ds_dy_from(y1000,x1000,y100,x100); %y1000....1st row of A..the observation of distance between 1000-100
    A(1,2) = ds_dx_from(y1000,x1000,y100,x100); %x1000
    A(1,7) = ds_dy_to(y1000,x1000,y100,x100); %y100
    A(1,8) = ds_dx_to(y1000,x1000,y100,x100); %x100
    
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
    
    A(10,5) = ds_dy_to(y100,x100,y3000,x3000);
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
    
    A(16,5) = ds_dy_to(y103,x103,y3000,x3000);
    A(16,6) = ds_dx_to(y103,x103,y3000,x3000);
    A(16,13) = ds_dy_from(y103,x103,y3000,x3000);
    A(16,14) = ds_dx_from(y103,x103,y3000,x3000);
    
    A(17,3) = ds_dy_to(y103,x103,y2000,x2000);
    A(17,4) = ds_dx_to(y103,x103,y2000,x2000);
    A(17,13) = ds_dy_from(y103,x103,y2000,x2000);
    A(17,14) = ds_dx_from(y103,x103,y2000,x2000);
    
    A(18,9) = ds_dy_to(y103,x103,y101,x101);
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
    
    A(36,9) = dr_dy_to(y103,x103,y101,x101);
    A(36,10) = dr_dx_to(y103,x103,y101,x101);
    A(36,13) = dr_dy_from(y103,x103,y101,x101);
    A(36,14) = dr_dx_from(y103,x103,y101,x101);
    A(36,21) = -1;
    
    A(37,1) = 1; %Partial derivatives of the equations for the coor of control points!!
    A(38,2) = 1; %1,2,3,4,5,6 because in these positions in the vector of unknowns are the control points
    A(39,3) = 1;
    A(40,4) = 1;
    A(41,5) = 1;
    A(42,6) = 1;
    
    %Normal matrix
    N = A'*P*A;
    rankA = rank(A);
    detN = det(N);
    %Vector of right hand side of normal equations
    n = A'*P*l;

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx = inv(N);
    
    %Solution of normal equations
    x_hat = Q_xx*n;
    
    %Adjusted unknowns
    X_hat = X_0+x_hat;
    
    %Update
    X_0 = X_hat;
    
    y1000=X_0(1);
    x1000=X_0(2);
    y2000=X_0(3);
    x2000=X_0(4);
    y3000=X_0(5);
    x3000=X_0(6);
    y100=X_0(7);
    x100=X_0(8);
    y101=X_0(9);
    x101=X_0(10);
    y102=X_0(11);
    x102=X_0(12);
    y103=X_0(13);
    x103=X_0(14);
    w1000=X_0(15);
    w2000=X_0(16);
    w3000=X_0(17);
    w100=X_0(18);
    w101=X_0(19);
    w102=X_0(20);
    w103=X_0(21);
    
    %Check 1
    max_x_hat = max(abs(x_hat)); %....of the solution of the normal equations...
	
	%Vector of residuals
    v = A*x_hat-l;
    v_gon = v(19:36,1)*200/pi;            %Convert to [gon]

    %Objective function
    vTPv = v'*P*v; 

    %Vector of adjusted observations
    L_hat = L+v;
    L_hat_gon = L_hat(19:36)*200/pi;       %Convert to [gon]

    % Distances
    Phi_X_hat(1) = dis(y1000,x1000,y100,x100);
    Phi_X_hat(2)=dis(y1000,x1000,y102,x102);
    Phi_X_hat(3)=dis(y2000,x2000,y103,x103);
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
    % Directions
    Phi_X_hat(19) = direction(y1000,x1000,y100,x100,w1000);
    Phi_X_hat(20)=direction(y1000,x1000,y102,x102,w1000);
    Phi_X_hat(21)=direction(y2000,x2000,y103,x103,w2000);
    Phi_X_hat(22)=direction(y2000,x2000,y101,x101,w2000);
    Phi_X_hat(23)=direction(y3000,x3000,y100,x100,w3000);
    Phi_X_hat(24)=direction(y3000,x3000,y103,x103,w3000);
    Phi_X_hat(25)=direction(y100,x100,y1000,x1000,w100);
    Phi_X_hat(26)=direction(y100,x100,y102,x102,w100);
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
    % Control points
    Phi_X_hat(37)= y1000;
    Phi_X_hat(38)= x1000;
    Phi_X_hat(39)= y2000;
    Phi_X_hat(40)= x2000;
    Phi_X_hat(41)= y3000;
    Phi_X_hat(42)= x3000;
    
    %Check 2
    Check2 = max(abs(L_hat-Phi_X_hat')); %final check:linearization of non-linear problem.
	
    %Update number of iterations
    iteration = iteration+1;

end

% Final Check
if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

%Convert to [gon] and check the quadrants
gon = X_0(15:21,1)*200/pi; %15-21 to the vector of unknowns
gon(1) = gon(1)+400;
gon(2) = gon(2)+400;
gon(3) = gon(3)+400;
gon(4) = gon(4)+400;
gon(5) = gon(5)+400;
gon(6) = gon(6)+400;
gon(7) = gon(7)+400;

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknows
s_X = sqrt(diag(S_XX_hat));
s_X_gon = s_X(15:21,1)*200/pi;           %Convert to [gon]

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));
s_L_hat_gon = s_L_hat(19:36)*200/pi; %Convert to [gon]....19-36 for directions

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));
s_v_gon = s_v(19:36,1)*200/pi;           %Convert to [gon]
 











