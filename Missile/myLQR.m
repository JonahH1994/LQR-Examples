function [t,S] = myLQR(A,B,Q,R,P1,tf )

% Create the time vector for ode45. This will do a backwards integration
tSpan =  [tf,0] ;

% Set the final state of the Riccatti equation to be P1
z0 = P1(:) ;
opts = odeset ;
sz = size(P1) ;

[t,x] = ode45( @dyn, tSpan, z0, opts, sz, A, B, Q, R ) ;
S = x ;

    function dz = dyn(~,z,sz,A,B,Q,R)
       
        % Since ode45 reshapes the matrix into a vector, reshape it into an
        % nxn matrix ;
        Pp = reshape(z(:),sz) ;
        
        % Calculate the derivative of the riccatti matrix:
        dP = -(Pp*A + A'*Pp -Pp*B*inv(R)*B'*Pp +Q) ;
        
        dz = dP(:) ;
        
    end

end