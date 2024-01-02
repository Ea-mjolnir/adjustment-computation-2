%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION
%        Gauss-Helmert model - Part II  
% 
%   Author         : Anastasia Pasioti
%   Version        : July 04, 2017
%   Last changes   : July 03, 2018
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Load file
plane = load('idun1.txt'); %we take all the coordinates

%Observations-MEASUREMENTS
x = plane(:,2);    %[m]
y = plane(:,3);    %[m]
z = plane(:,4);    %[m]

%Vector of observations
l = [x; y; z]; 

%Number of observations / number of POINTS
no_n = length(l)/3 ;  %4 observations

%Initial values for the unknowns
nx = 0.5774; %mono OXI 0...the same result for values=3,10,50!!
ny = 0.5774;
nz = 0.5774;
d = 0;

%Vector of initial values for the unknowns
X_0 = [nx ny nz d]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u+1;  %+1 BECAUSE WE HAVE ONE CONTRAINT EQUATION!!! 

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
S_ll = eye(3*no_n); %3 because the coordinates of points (X, Y, Z)!!

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_ll = 1/sigma_0^2*S_ll;

%Weight matrix
P = inv(Q_ll);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-10; %FOR THE COMPUTATION ERROR
delta = 10^-12; %FOR LINEARIZATION ERROR
max_psi = 10^Inf; %MAXIMUM NUMBER OF Ø:CONDITION EQUATION

%Initialization
psi = zeros(no_n,1); %4 CONDITION EQUATIONS = 4 POINTS
A = zeros(no_n,no_u); %4x4
B = zeros(no_n,3*no_n); %4x12

%Initial values for the residuals
v = zeros(3*no_n,1); %we give the value zero as vector!!!% 3 is x y z

%Number of iterations
iteration = 0;

%Constants
c1 = zeros(no_n,1);
c2 = 1; %from the constrain equations!!   ||n||=1

while max_psi>epsilon || Check2>delta                
    
    %Psi function...condition equation of the Adjustment
      psi = nx*(x+v(1:end/3))+ny*(y+v(end/3+1:end/1.5))+nz*(z+v(end/1.5+1:end))-d;

    %Gamma function
      gamma = nx^2+ny^2+nz^2;
    
     %Matrices with the elements from the Jacobian matrices J1, J2 and J3
      A = [x+v(1:end/3) y+v(end/3+1:end/1.5) z+v(end/1.5+1,end) -ones(no_n,1)];
    
      B = [diag(nx*ones(no_n,1)) diag(ny*ones(no_n,1)) diag(nz*ones(no_n,1))];

      C = [2*nx 2*ny 2*nz 0]; %from the Constraint equation!!!!

    %Vectors of misclosures...........................................................
     w1 = -B*v+psi-c1; %zeros(no_n,1)
     w2 = gamma-c2;
    
    %Normal matrix
     N_ext = [-A'*inv(B*Q_ll*B')*A C';C 0];
     
    %Vector of right hand side of normal equations
     n_ext = [A'*inv(B*Q_ll*B')*w1; -w2];
    
    %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx_ext = inv(N_ext);
    
    %Solution of the normal equations
     x_hat = Q_xx_ext*n_ext;
    
    %Langrange multipliers
     k1 = inv(B*Q_ll*B')*(-A*x_hat(1:end-1)-w1); %the last entry of x_hat is the vector of k2
     k2 = x_hat(end);
    
    %Vector of residuals
     v_new = Q_ll*B'*k1;
       
    %Update
     X_0 =  X_0+x_hat(1:no_u); %Because we do not want the k %Here Adj.model with conditions and constraints the vector of unknows is at the top!
    
    nx = X_0(1);
    ny = X_0(2);
    nz = X_0(3);
    d = X_0(4);

    %Update residuals
     v = v_new;
    
    %Check 1
    max_psi = max(abs(psi-c1));
    
    %vTPv function
     vTPv = v'*P*v;
 
    %Vector of adjusted observations
     l_hat = l+v;

    %Check 2
     Check2 = max(abs(nx*l_hat(1:end/3)+ny*l_hat(end/3+1:end/1.5)+nz*l_hat(end/1.5+1:end)-d));                                 
     %We need the max of the absolute values of the condition equations
    
    %Update number of iterations
    iteration = iteration+1;
    
     
end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);    %a posteriori 

%VC matrix of adjusted unknowns
S_XX_hat = -s_0^2*Q_xx_ext(1:no_u,1:no_u);

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));

%Cofactor matrix of the residuals
Q_vv = Q_ll*B'*inv(B*Q_ll*B')*B*Q_ll;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));

%Cofactor matrix of adjusted observations
Q_ll_hat = Q_ll-Q_vv;

%VC matrix of adjusted observations
S_ll_hat = s_0^2*Q_ll_hat;

%Standard deviation of the adjusted observations
s_l_hat = sqrt(diag(S_ll_hat));

%Results for the unknowns
nx 
ny 
nz 
d 
%%
%Load file
plane = load('idun1.txt'); %we take all the coordinates
nx=0.0037;		
ny=0.0057;	
nz=0.99997;	
d=-99.8258;

x = plane(:,2);  %[m]
y = plane(:,3);  %[m]
z = plane(:,4);  %[m]

figure
%Plot points
plot3(x,y,z,'or')
hold on

%Generate x and y data
[x1,y1] = meshgrid(95:1:110); %5 is the step

%Solve for z data
%z1 = (-1*(nx.*x1+ny.*y1+d.*ones(8,1)))/nz;
z1 = -1/nz*(nx.*x1+ny.*y1+d);

%Plot the plane/surface
surf(x1,y1,z1)

xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Observed points and adjusted plane.')

% Data to write to the file
X_new=reshape(x1,[],1);
Y_new=reshape(y1,[],1);
Z_new=reshape(z1,[],1);

coord_plain=[X_new,Y_new,Z_new];
%fileID = fopen('coord_plain.txt','w');
%fprintf(fileID,'%4d %4d %4d\n',coord_plain);
%fclose(fileID);

%%
x2 =coord_plain(:,1);  %[m]
y2 =coord_plain(:,2);  %[m]
z2 =coord_plain(:,3);  %[m]
%Plot the plane/surface
plot3(x2,y2,z2,'r')





















