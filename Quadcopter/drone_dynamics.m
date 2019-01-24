%% Deriver for the drone dynamics:

clear all ;

syms m g Ixx Iyy Izz H S t ;
% S = psi star or reference heading

x = sym('x',[12,1]) ;
u = sym('u',[4,1]) ;

eq1  = x(2) ;
eq2  = (sin(x(11))*sin(x(7)) + cos(x(11))*sin(x(9))*cos(x(7)))*u(1)/m ;
eq3  = x(4) ;
eq4  = (-cos(x(11))*sin(x(7)) + sin(x(11))*sin(x(9))*cos(x(7)))*u(1)/m ;
eq5  = x(6) ;
eq6  = -g + cos(x(9))*cos(x(7))*u(1)/m ;
eq7  = x(8) ;
eq8  = (Iyy-Izz)/Ixx*x(10)*x(8) + u(2)/Ixx ;
eq9  = x(10) ;
eq10 = (Izz-Ixx)/Iyy*x(8)*x(12) + u(3)/Iyy ;
eq11 = x(12) ;
eq12 = (Ixx-Iyy)/Izz*x(8)*x(10) + u(4)/Izz ;

eqns = [eq1; eq2; eq3; eq4; eq5; eq6; eq7; eq8; eq9; eq10; eq11; eq12 ] ;

x0 = [0;0;0;0;H;0;0;0;0;0;S;0] ;
u0 = [m*g;0;0;0] ;

uA = jacobian( eqns, x ) ;
uB = jacobian( eqns, u ) ;
A = subs(uA, [x,u], [x0,u0] ) ;
B = subs(uB, [x,u], [x0,u0] ) ;

matlabFunction( A, 'File', 'Amatrix' ) ;
matlabFunction( B, 'File', 'Bmatrix' ) ;

% transA = expm(A*t) ;
% transAt = expm(A.'*t) ;
% 
% integralArg = transAt * transA ;