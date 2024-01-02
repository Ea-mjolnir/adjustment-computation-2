%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION
%        Gauss-Helmert model - Part II  
% 
%   Author         : Elina
%   Version        : July 04, 2017
%   Last changes   : July 03, 2018
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format long g;

%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------

%Observations
X1 = 4540134.2780; %[m] X-coor of TARGET-GLOBAL SYSTEM!!
Y1 = 382379.8964; %[m] Y-coor of TARGET-GLOBAL SYSTEM!!
x1 = 4540124.0940; %[m] x-coor of START-LOCAL SYSTEM!!
y1 = 382385.9980; %[m] y-coor of START-LOCAL SYSTEM!!

X2 = 4539937.3890;
Y2 = 382629.7872;
x2 = 4539927.2250;
y2 = 382635.8691;

X3 = 4539979.7390;
Y3 = 381951.4785;
x3 = 4539969.5670;
y3 = 381957.5705;

X4 = 4540326.4610;
Y4 = 381895.0089;
x4 = 4540316.2940;
y4 = 381901.0932;

X5 = 4539216.3870;
Y5 = 382184.4352;
x5 = 4539206.2110;
y5 = 382190.5278;

%Vector of observations
l = [X1 Y1 X2 Y2 X3 Y3 X4 Y4 X5 Y5 x1 y1 x2 y2 x3 y3 x4 y4 x5 y5]';

%Centroids
X_c = (X1+X2+X3+X4+X5)/5;
Y_c = (Y1+Y2+Y3+Y4+Y5)/5;
x_c = (x1+x2+x3+x4+x5)/5;
y_c = (y1+y2+y3+y4+y5)/5;

%Coordinates reduced to the centroid
X1 = X1 - X_c;
Y1 = Y1 - Y_c;
X2 = X2 - X_c;
Y2 = Y2 - Y_c;
X3 = X3 - X_c;
Y3 = Y3 - Y_c;
X4 = X4 - X_c;
Y4 = Y4 - Y_c;
X5 = X5 - X_c;
Y5 = Y5 - Y_c;

x1 = x1 - x_c;
y1 = y1 - y_c;
x2 = x2 - x_c;
y2 = y2 - y_c;
x3 = x3 - x_c;
y3 = y3 - y_c;
x4 = x4 - x_c;
y4 = y4 - y_c;
x5 = x5 - x_c;
y5 = y5 - y_c;

%Vector of observations.
l = [X1 Y1 X2 Y2 X3 Y3 X4 Y4 X5 Y5 x1 y1 x2 y2 x3 y3 x4 y4 x5 y5]';

%Number of observations / points....10 (5 points in one system + 5 points in the other system)
no_n = length(l)/2; %Not all the observations!!! THE NUMBER OF POINTS - NUMBER OF CONDITION EQUATIONS!!!

%Initial values for the unknowns
x0 = X_c - x_c;
y0 = Y_c - y_c;
m = 0.8999;
alpha = 0.002;   %[rad]

%for linearization.
a = m*cos(alpha);
o = m*sin(alpha);

%Vector of initial values for the unknowns.
X_0 = [a o x0 y0]'; % a because of linearization

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u;   %withoout constraint equations

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%Weight matrix...%Firstly the 10 values of P_XY AND follow the ten values of p_xy
P_XY = diag([10.0000 14.2857 0.8929 1.4286 7.1429 10.0000 2.2222 3.2259 7.6923 11.1111]);
P_xy = diag([5.8824 12.5000 0.9009 1.7241 7.6923 16.6667 4.1667 6.6667 8.3333 16.6667]);
P = [P_XY zeros(no_n, no_n)
     zeros(no_n, no_n) P_xy];
    
Q_ll = inv(P);    

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-10; %FOR THE COMPUTATION ERROR
delta = 10^-12; %FOR LINEARIZATION ERROR
max_psi = 10^Inf;

%Initialization
psi = zeros(no_n,1); %10 CONDITION EQUATIONS (5 points * 2)
A = zeros(no_n,no_u); %10x4...wrt to the unknowns
B = zeros(no_n,2*no_n); %10x20....wrt to the residuals

%Initial values for the residuals
v = zeros(2*no_n,1); %we give the value zero as vector!!!

%Number of iterations
iteration = 0;

%Constants
c1 = zeros(no_n,1);

while max_psi>epsilon || Check2>delta
    
    % (X+vx) -a*(x+vx)+o*(y+vy)-x_o;
    % (Y+vy) -o*(x+vx)-a*(y+vy)-y_o;
    
    %Psi function
     psi = [(l(1)+v(1))-a*(l(11)+v(11))+o*(l(12)+v(12))-x0;
            (l(2)+v(2))-o*(l(11)+v(11))-a*(l(12)+v(12))-y0;
            
            (l(3)+v(3))-a*(l(13)+v(13))+o*(l(14)+v(14))-x0;
            (l(4)+v(4))-o*(l(13)+v(13))-a*(l(14)+v(14))-y0;
            
    % (X+vx) -a*(x+vx)+o*(y+vy)-x_o;
    % (Y+vy) -o*(x+vx)-a*(y+vy)-y_o;            
            
            (l(5)+v(5))-a*(l(15)+v(15))+o*(l(16)+v(16))-x0;
            (l(6)+v(6))-o*(l(15)+v(15))-a*(l(16)+v(16))-y0;
            
            (l(7)+v(7))-a*(l(17)+v(17))+o*(l(18)+v(18))-x0;
            (l(8)+v(8))-o*(l(17)+v(17))-a*(l(18)+v(18))-y0;
           
            (l(9)+v(9))-a*(l(19)+v(19))+o*(l(20)+v(20))-x0;
            (l(10)+v(10))-o*(l(19)+v(19))-a*(l(20)+v(20))-y0];
    
     %Matrices with the elements from the Jacobian matrices J1, J2.
     % (X+vx) -a*(x+vx)+o*(y+vy)-x_o;
     % (Y+vy) -o*(x+vx)-a*(y+vy)-y_o; 
     
     %           a             o       x0 y0
     A = [-(l(11)+v(11)) (l(12)+v(12)) -1 0;
          -(l(12)+v(12)) -(l(11)+v(11)) 0 -1;
          
          -(l(13)+v(13)) (l(14)+v(14)) -1 0;
          -(l(14)+v(14)) -(l(13)+v(13)) 0 -1;
          
          -(l(15)+v(15)) (l(16)+v(16)) -1 0;
          -(l(16)+v(16)) -(l(15)+v(15)) 0 -1;
          
          -(l(17)+v(17)) (l(18)+v(18)) -1 0;
          -(l(18)+v(18)) -(l(17)+v(17)) 0 -1;
          
          -(l(19)+v(19)) (l(20)+v(20)) -1 0;
          -(l(20)+v(20)) -(l(19)+v(19)) 0 -1];
      
    % (X+vx) -a*(x+vx)+o*(y+vy)-x_o;
    % (Y+vy) -o*(x+vx)-a*(y+vy)-y_o;       
    
     B = [diag(ones(no_n,1))];
     B(1,11) = -a;
     B(1,12) = o;
     B(2,11) = -o;
     B(2,12) = -a;
     
     B(3,13) = -a;
     B(3,14) = o;
     B(4,13) = -o;
     B(4,14) = -a;
     
     B(5,15) = -a;
     B(5,16) = o;
     B(6,15) = -o;
     B(6,16) = -a;
     
     B(7,17) = -a;
     B(7,18) = o;
     B(8,17) = -o;
     B(8,18) = -a;
     
     B(9,19) = -a;
     B(9,20) = o;
     B(10,19) = -o;
     B(10,20) = -a;
     
    %Vectors of misclosures.
    w1 = -B*v+psi-c1; %C=0
    
    %Normal matrix.
    N_ext = [B*Q_ll*B'     A; 
                 A' zeros(no_u,no_u)]; 
             
    %like in exercise 9 where we do NOT have constraints.
     
    %Vector of right hand side of normal equations.
    n_ext = [-w1; zeros(no_u,1)];
    
    %Inversion of normal matrix / Cofactor matrix of the unknowns.
    Q_xx_ext = inv(N_ext);
    
    %Solution of the normal equations.
    x_hat = Q_xx_ext*n_ext;
    
    %Langrange multipliers.
    k = inv(B*Q_ll*B')*(-w1);
    
    %Vector of residuals......................until here the residuals have the value 0!!
    v_new = Q_ll*B'*k; %or k or x_hat(1:no_n)
       
    %Update.
    X_0 = X_0 + x_hat(end-no_u+1:end); %We take the last 10 values of vector [k X]...
    
    a = X_0(1);
    o = X_0(2);
    x0 = X_0(3);
    y0 = X_0(4);

    %Update residuals
    v = v_new;
    
    %Check 1
    max_psi = max(abs(psi));
    
    %vTPv function...aposteriori std of the Adjustement
    vTPv = v'*P*v;
 
    %Vector of adjusted observations
    l_hat = l+v;   % residuals have been updated v=v_new;
    
    psi = [l_hat(1)-a*l_hat(11)+o*l_hat(12)-x0; %I create again the condition equations by the adjusted observations!!!
           l_hat(2)-o*l_hat(11)-a*l_hat(12)-y0;
            
            l_hat(3)-a*l_hat(13)+o*l_hat(14)-x0;
            l_hat(4)-o*l_hat(13)-a*l_hat(14)-y0;
           
            l_hat(5)-a*l_hat(15)+o*l_hat(16)-x0;
            l_hat(6)-o*l_hat(15)-a*l_hat(16)-y0;
            
            l_hat(7)-a*l_hat(17)+o*l_hat(18)-x0;
            l_hat(8)-o*l_hat(17)-a*l_hat(18)-y0;
           
            l_hat(9)-a*l_hat(19)+o*l_hat(20)-x0;
            l_hat(10)-o*l_hat(19)-a*l_hat(20)-y0];
    
    %Check 2
    Check2 = max(abs(psi));
    
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
%S_XX_hat = -s_0^2*Q_xx_ext(1:no_u,1:no_u);
S_XX_hat = -s_0^2*Q_xx_ext(end-no_u+1:end, end-no_u+1:end);

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));

%Cofactor matrix of the residuals
%Q_vv = Q_ll*B'*inv(B*Q_ll*B')*B*Q_ll;
Q_vv = Q_ll*B'*Q_xx_ext(1:no_n,1:no_n)*B*Q_ll;

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

%De-reduced unknowns
x0 = X_0(3)+X_c;
y0 = X_0(4)+Y_c;

%De-reduced coordinates
X1 = l_hat(1)+X_c;
Y1 = l_hat(2)+Y_c;
X2 = l_hat(3)+X_c;
Y2 = l_hat(4)+Y_c;
X3 = l_hat(5)+X_c;
Y3 = l_hat(6)+Y_c;
X4 = l_hat(7)+X_c;
Y4 = l_hat(8)+Y_c;
X5 = l_hat(9)+X_c;
Y5 = l_hat(10)+y_c;
x1 = l_hat(11)+x_c;
y1 = l_hat(12)+y_c;
x2 = l_hat(13)+x_c;
y2 = l_hat(14)+y_c;
x3 = l_hat(15)+x_c;
y3 = l_hat(16)+y_c;
x4 = l_hat(17)+x_c;
y4 = l_hat(18)+y_c;
x5 = l_hat(19)+x_c;
y5 = l_hat(20)+y_c;

l_final = [X1 Y1 X2 Y2 X3 Y3 X4 Y4 X5 Y5 x1 y1 x2 y2 x3 y3 x4 y4 x5 y5]';

%Results for the unknowns.
a;
o;
x0;
y0;

%compute m and alpha.
alpha = atan(o/a)*200/pi;
m = a/cos(alpha);










