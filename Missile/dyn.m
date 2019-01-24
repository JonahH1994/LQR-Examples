function dz = dyn(t,z,A,B,R,S,T)

%K = interp( T, S, t ) ;
[~,b1] = min(abs(T-t));
K = R\B'*reshape(S(b1,:),size(A)) ;
Aa = A - B*K ;

dz = Aa * z ;

end