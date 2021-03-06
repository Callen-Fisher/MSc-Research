function predict = predictFunctionFullEKF(Tp,Ty,dph,dth,ph,signdph,signdth,th)
%PREDICTFUNCTIONFULLEKF
%    PREDICT = PREDICTFUNCTIONFULLEKF(TP,TY,DPH,DTH,PH,SIGNDPH,SIGNDTH,TH)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    02-Jul-2015 12:38:37

t2 = cos(ph);
t3 = t2.^2;
t4 = cos(th);
t5 = dph.^2;
t6 = dth.^2;
t7 = ph.*2.0;
t8 = sin(t7);
predict = [0.0,0.0,(Tp.*6.666666666666667e4+t4.*6.468e4-t5.*sin(th.*2.0).*5.445e2+dph.*dth.*t8.*1.52e3-signdth.*t3.*t6.*2.90550645e2)./(t3.*1.52e3+1.089e3),-(Ty.*(-6.666666666666667e4)+t6.*t8.*7.6e2+t4.*(signdph.*t5.*2.90550645e2-dph.*dth.*sin(th).*2.178e3))./(t4.^2.*1.089e3+1.52e3),dth,dph];
