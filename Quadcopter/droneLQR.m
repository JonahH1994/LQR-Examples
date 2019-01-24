function droneLQR()

% Define constants:
g = 9.81 ; % gravity 
m = 2.1 ; % mass of quadcopter
Ixx = 2.1385 ; % Inertia about x axis
Iyy = Ixx ; % Inertia bout y axis
Izz = 3.7479 ; % Inertia about z axis
tf = 8 ; % Final time

% Motor specific constants were taken from the quadrotor constants matlab
% file provided for the final project.

%% First set of initial conditions
% Define desired altitude H* and heading angle psi*:
H = -30 ;
psid = 4*pi/3 ;
x0 = [0;0;0;0;H;0;0;0;0;0;psid;0] ;

% Define A and B matrices:
A = Amatrix(psid,g) ;
B = Bmatrix(Ixx,Iyy,Izz,m) ;

% Define Q and R matrices for LQR:
weight = 10^6 ;
Q = eye(size(A)) * weight ;
R = eye(size(B,2),size(B,2)) * 1/weight ;

% Find the weighting matrix for inputs using lqr:
[K,~,~] = lqr( A, B, Q, R ) ;

% Define closed loop state:
Ad = A-B*K ;
Bd = B ;
% output the error in the altitude, heading and angles phi and theta:
C = [0,0,0,0,1,0,0,0,0,0,0,0; ...
     0,0,0,0,0,0,0,0,0,0,1,0; ...
     0,0,0,0,0,0,1,0,0,0,0,0; ...
     0,0,0,0,0,0,0,0,1,0,0,0 ] ;
D = 0 ;
sysCL = ss(Ad,Bd,C,D) ;

% Simulate step response:
tin = 0:0.1:tf;
U = zeros(size(tin,2),size(B,2)) ;
[y,t,~] = lsim( sysCL,U,tin,x0) ;

figure
plot( t, y(:,1) ) ;
hold on ;
plot( t, y(:,2) ) ;
legend( 'Altitude', 'Heading' ) ;
xlabel( 'Time (s)' )  ;
ylabel( 'Error' ) ;
title( 'H^*=-30 & \psi^* = 4\pi/3' ) ;

%% Second set of initial conditions:

% Define desired altitude H* and heading angle psi*:
H = 20 ;
psid = -5*pi/6 ;
x0 = [0;0;0;0;H;0;0;0;0;0;psid;0] ;

% Define A and B matrices:
A = Amatrix(psid,g) ;
B = Bmatrix(Ixx,Iyy,Izz,m) ;

% Run lqr:
[K,~,~] = lqr( A, B, Q, R ) ;

% Define closed loop system:
Ad = A - B*K ;
sysCL = ss(Ad,Bd,C,D) ;

% Output:
[y,t,~] = lsim( sysCL,U,tin,x0) ;

figure
plot( t, y(:,1) ) ;
hold on ;
plot( t, y(:,2) ) ;
legend( 'Altitude', 'Heading' ) ;
xlabel( 'Time (s)' )  ;
ylabel( 'Error' ) ;
title( 'H^*=20 & \psi^* = -5\pi/6' ) ;

%% Third set of initial conditions:

% Define desired altitude H* and heading angle psi*:
H = -30 ;
psid = 4*pi/3 ;
x0 = [0;0;0;0;H;0;0;0;0;0;psid;0] ;

% Define A and B matrices:
A = Amatrix(psid,g) ;
B = Bmatrix(Ixx,Iyy,Izz,m) ;

% Run lqr:
[K,~,~] = lqr( A, B, Q, R ) ;

% Define closed loop system:
Ad = A - B*K ;
sysCL = ss(Ad,Bd,C,D) ;

% Output:
U = zeros(size(tin,2),size(B,2)) ;
% Impulse disturbances
imp = 100000 ;
u1 = tin == 3 ;
U(:,1) = u1' * imp ;
u2 = tin == 4;
U(:,2) = u2' * imp ;
u3 = tin == 2 ;
U(:,3) = u3' * imp;
u4 = tin == 3.5 ;
U(:,4) = u4' *imp ;
[y,t,~] = lsim( sysCL,U,tin,x0) ;

figure
plot( t, y(:,1) ) ;
hold on ;
plot( t, y(:,2) ) ;
legend( 'Altitude', 'Heading' ) ;
xlabel( 'Time (s)' )  ;
ylabel( 'Error' ) ;
title( 'H^*=-30, \psi^* = 4\pi/3 & impulses in U' ) ;

% The impulses had no effect on the output of the desired values but
% impacted the values that we had not specified, however, the error still 
% remained fairly small:

figure
plot( t, y(:,3) ) ;
hold on ;
plot( t, y(:,4) ) ;
legend( '\phi', '\theta' ) ;
xlabel( 'Time (s)' ) ;
ylabel( 'Error' ) ;
title( 'H^*=-30, \psi^* = 4\pi/3 & impulses in U' ) ;

%% Fourth set of initial conditions:

% Define desired altitude H* and heading angle psi*:
H = 20 ;
psid = -5*pi/6 ;
x0 = [0;0;0;0;H;0;0;0;0;0;psid;0] ;

% Define A and B matrices:
A = Amatrix(psid,g) ;
B = Bmatrix(Ixx,Iyy,Izz,m) ;

% Run lqr:
[K,~,~] = lqr( A, B, Q, R ) ;

% Define closed loop system:
Ad = A - B*K ;
sysCL = ss(Ad,Bd,C,D) ;

% Output:
U = zeros(size(tin,2),size(B,2)) ;
% Impulse disturbances
imp = 100000 ;
u1 = tin == 3 ;
U(:,1) = u1' * imp ;
u2 = tin == 4;
U(:,2) = u2' * imp ;
u3 = tin == 2 ;
U(:,3) = u3' * imp;
u4 = tin == 3.5 ;
U(:,4) = u4' *imp ;
[y,t,~] = lsim( sysCL,U,tin,x0) ;

figure
plot( t, y(:,1) ) ;
hold on ;
plot( t, y(:,2) ) ;
legend( 'Altitude', 'Heading' ) ;
xlabel( 'Time (s)' )  ;
ylabel( 'Error' ) ;
title( 'H^*=20, \psi^* = -5\pi/6 & impulses in U' ) ;

end