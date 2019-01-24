clear

tf = 500 ;

%Missile Parameters:
S = 0.44;
rho = 973.3;
a = 1036.4;
m = 13.98;
g = 32.2;
l = 0.75;
Iy = 182.5;

%Enter Mach Number (2-4):
M = 2;
V = M*a;
Q = 0.7*rho*(M^2);

%Missile Aerodynamic Coefficients:
n1 = -0.1696*(1 + 1/3*(3-M));
n0 = -0.34;
m1 = 0.051*(1 + 8/3*(M-3));
m0 = -0.206;

%Missile Stability Coefficients:
Zalpha = Q*S*n0/m;
Zdelta = Q*S*n1/m;
Malpha = Q*S*l*m1/Iy;
Mdelta = Q*S*l*m0/Iy;

%Missile State-Space Matrices:
A = [Zalpha/V 1; Malpha 0];
B = [Zdelta/V; Mdelta];
C = [Zalpha 0; 0 1];
D = [Zdelta; 0];

% To get the natural frequency of the target, multiply the wave number by
% the velocity:
wn = 2 * pi/5280 * V ;

% Initial y position of the missile
y0 = -20000 ;

format long ;

% The state matrices of my system:
An = [   0, 1,           0, 0,     0, 0 ; ...
         0, 0, Q*S*l*m1/Iy, 0,     0, 0 ; ...
         0, 1,  Q*S*n0/m/V, 0,     0, 0 ; ...
         V, 0,          -V, 0,     0,-1 ; ...
         0, 0,           0, 0,     0, 1 ; ...
         0, 0,           0, 0, -wn^2, 0] ;

Bn = [ 0; Q*S*l*m0/Iy; Q*S*n1/m/V; 0; 0; 0] ;

% I would like to output y, and the target trajectory of y
Cn = [0,0,0,1, 1, 0 ;0,0,0,0,1,0;zeros(1,6)] ;
Dn = [0;0;0] ;

% Define the weights for the LQR:
R = 1e10 ;
Q = [zeros(1,6); zeros(1,6); zeros(1,6); ...
   0, 0, 0, 1, 0, 0; zeros(1,6); zeros(1,6) ] ; 
Q = Q * 2e1 ;

% Call the LQR function which will return the riccatti matrix at each time
% step from t = 0 to t = tf:
[t11,S1] = myLQR( An,Bn,Q,R,zeros(size(An)),tf ) ;

% Define the initial conditions for the simulation
t = 0;
r = 1000*(-1/5*sin(2*pi*(t-10)/5280)) ;
dr = -1/5*1000*2*pi/5280*cos(10*2*pi/5280)*V ;
x0 = [0;0;0;y0-r(1);r(1);dr] ;

% ode45 is run here using RHS (Right Hand Side) file dyn which computes the
% state of the sytem at each time step. Note that the dyn function also
% takes as inputs the A and B state matrices, the weight R and ricatti
% matrices at each time stem to calculate K as well as the time vector that
% corresponds to the ricatti matrix.
opts = odeset ;
tSpan = [0 tf] ;
[t1,xx] = ode45( @dyn, tSpan, x0, opts, An, Bn,R, S1, t11' ) ;

% Extract the desired outputs:
% I use the final K value to output the control input into the system
% since u = -Kx: (I could be more robust and keep delta in the state but
% I determined that is is reasonable to use this K value). This is the
% value that lqr in matlab uses.
K = R\Bn'*reshape(S1(end,:),size(An)) ;
Cc = [Cn(1:2,:);-K];
y1 = (Cc * xx')' ;

% Plot the trajectory of the target and of the missile:
figure
plot( t1*V, y1(:,1), 'k', 'LineWidth', 2 ) ;
hold on ;
plot( t1*V, y1(:,2), '--k' ) ;
xlabel( 'X (m)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Missile and Target Trajectory (Base Case)' ) ;
legend( 'Missile', 'Target' ) ;
hold off ;

figure
plot( t1*V, ([0,0,0,1,0,0]*xx')' ) ;
hold on ;
xlabel( 'X (m)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Error in the Missiles Trajectory (Base Case)' ) ;
hold off ;

figure
plot( t1*V, y1(:,3) ) ;
hold on ;
xlabel( 'X (m)' ) ;
ylabel( '\Delta (rad)' ) ;
title( 'Plot of the Missiles Delta Angle Along the Trajectory' ) ;
hold off ;

figure
plot( t1*V, ([1,0,0,0,0,0]*xx')' ) ;
hold on ;
xlabel( 'X (m)' ) ;
ylabel( '\theta (rad)' ) ;
title( 'Plot of the Missiles Orientation (\theta) Along the Trajectory (Base Case)' ) ;
hold off ;

figure
plot( t1, y1(:,1) ) ;
hold on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Missiles Trajectory v. Time (Base Case)' ) ;
hold off ;

%% Case Studies:

%% Case Study 1:

% Start the missle at an orientation of 5 degrees:
x0 = [10*pi/180;0;0;y0-r(1);r(1);dr] ;

% Run ode45 and extract outputs:
[t1,xx] = ode45( @dyn, tSpan, x0, opts, An, Bn,R, S1, t11' ) ;
Cc = Cn(1:2,:);
y1 = (Cc * xx')' ;

% Plot the trajectory of the target and of the missile:
figure
plot( t1*V, y1(:,1), 'k', 'LineWidth', 2 ) ;
hold on ;
plot( t1*V, y1(:,2), '--k' ) ;
xlabel( 'X (m)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Missile and Target Trajectory (Case Study 1)' ) ;
legend( 'Missile', 'Target' ) ;
hold off ;

figure
plot( t1*V, ([0,0,0,1,0,0]*xx')' ) ;
hold on ;
xlabel( 'X (m)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Error in the Missiles Trajectory (Case Study 1)' ) ;
hold off ;

figure
plot( t1*V, ([1,0,0,0,0,0]*xx')' ) ;
hold on ;
xlabel( 'X (m)' ) ;
ylabel( '\theta (rad)' ) ;
title( 'Plot of the Missiles Orientation (\theta) Along the Trajectory (Case Study 1)' ) ;
hold off ;

figure
plot( t1, y1(:,1) ) ;
hold on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Missiles Trajectory v. Time (Case Study 1)' ) ;
hold off ;

%% Case Study 1:
% Start the missle at a y position of 20,000:
y0 = 30000 ;
x0 = [0;0;0;y0-r(1);r(1);dr] ;

% Run ode45 and extract outputs:
[t1,xx] = ode45( @dyn, tSpan, x0, opts, An, Bn,R, S1, t11' ) ;
Cc = Cn(1:2,:);
y1 = (Cc * xx')' ;

[t1,xx] = ode45( @dyn, tSpan, x0, opts, An, Bn,R, S1, t11' ) ;

Cc = Cn(1:2,:);
y1 = (Cc * xx')' ;

% Plot the trajectory of the target and of the missile:
figure
plot( t1*V, y1(:,1), 'k', 'LineWidth', 2 ) ;
hold on ;
plot( t1*V, y1(:,2), '--k' ) ;
xlabel( 'X (m)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Missile and Target Trajectory (Case Study 2)' ) ;
legend( 'Missile', 'Target' ) ;
hold off ;

figure
plot( t1*V, ([0,0,0,1,0,0]*xx')' ) ;
hold on ;
xlabel( 'X (m)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Error in the Missiles Trajectory (Case Study 2)' ) ;
hold off ;

figure
plot( t1*V, ([1,0,0,0,0,0]*xx')' ) ;
hold on ;
xlabel( 'X (m)' ) ;
ylabel( '\theta (rad)' ) ;
title( 'Plot of the Missiles Orientation (\theta) Along the Trajectory (Case Study 2)' ) ;
hold off ;

figure
plot( t1, y1(:,1) ) ;
hold on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Missiles Trajectory v. Time (Case Study 2)' ) ;
hold off ;

%% Case Study 3:
V = 2 * V ;

% Need to redefine the state:
S = 0.44;
Q = 0.7*rho*(M^2);
wn = 2 * pi/5280 * V ;

% The state matrices of my system:
An = [   0, 1,           0, 0,     0, 0 ; ...
         0, 0, Q*S*l*m1/Iy, 0,     0, 0 ; ...
         0, 1,  Q*S*n0/m/V, 0,     0, 0 ; ...
         V, 0,          -V, 0,     0,-1 ; ...
         0, 0,           0, 0,     0, 1 ; ...
         0, 0,           0, 0, -wn^2, 0] ;

Bn = [ 0; Q*S*l*m0/Iy; Q*S*n1/m/V; 0; 0; 0] ;

Q = [zeros(1,6); zeros(1,6); zeros(1,6); ...
   0, 0, 0, 1, 0, 0; zeros(1,6); zeros(1,6) ] ; 
Q = Q * 2e1 ;

[t11,S1] = myLQR( An,Bn,Q,R,zeros(size(An)),tf ) ;

dr = -1/5*1000*2*pi/5280*cos(10*2*pi/5280)*V ;
x0 = [0;0;0;y0-r(1);r(1);dr] ;

[t1,xx] = ode45( @dyn, tSpan, x0, opts, An, Bn,R, S1, t11' ) ;

Cc = Cn(1:2,:);
y1 = (Cc * xx')' ;

% Plot the trajectory of the target and of the missile:
figure
plot( t1*V, y1(:,1), 'k', 'LineWidth', 2 ) ;
hold on ;
plot( t1*V, y1(:,2), '--k' ) ;
xlabel( 'X (m)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Missile and Target Trajectory (Case Study 3)' ) ;
legend( 'Missile', 'Target' ) ;
hold off ;

figure
plot( t1*V, ([0,0,0,1,0,0]*xx')' ) ;
hold on ;
xlabel( 'X (m)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Error in the Missiles Trajectory (Case Study 3)' ) ;
hold off ;

figure
plot( t1*V, ([1,0,0,0,0,0]*xx')' ) ;
hold on ;
xlabel( 'X (m)' ) ;
ylabel( '\theta (rad)' ) ;
title( 'Plot of the Missiles Orientation (\theta) Along the Trajectory (Case Study 3)' ) ;
hold off ;

figure
plot( t1, y1(:,1) ) ;
hold on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Y (m)' ) ;
title( 'Plot of the Missiles Trajectory v. Time (Case Study 3)' ) ;
hold off ;