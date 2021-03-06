function predict = predictFunctionFullEKF(Tp1,Tp2,dth1,dth2,signdth1,signdth2,th1,th2)
%PREDICTFUNCTIONFULLEKF
%    PREDICT = PREDICTFUNCTIONFULLEKF(TP1,TP2,DTH1,DTH2,SIGNDTH1,SIGNDTH2,TH1,TH2)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    23-May-2015 11:38:31

t2 = dth1.^2;
t3 = th1.*2.0;
t8 = th2.*2.0;
t4 = t3-t8;
t5 = th1-th2;
t6 = dth2.^2;
t7 = cos(t5);
t9 = cos(t4);
t10 = t9.*1.458e3;
t11 = t10-2.1985e4;
t12 = 1.0./t11;
t13 = sin(t5);
t14 = sin(t4);
predict = [0.0,0.0,t12.*(Tp1.*(-6.566666666666667e5)+cos(th1-th2.*2.0).*3.969e4-cos(th1).*4.81572e5+Tp2.*t7.*3.6e5+signdth1.*t2.*2.53333332e2+t2.*t14.*1.458e3+t6.*t13.*5.319e3+signdth2.*t6.*t7.*1.28496375e2),-t12.*(Tp2.*1.586666666666667e6-cos(t3-th2).*1.42884e5+cos(th2).*2.06976e5-Tp1.*t7.*3.6e5+signdth2.*t6.*5.66335875e2+t2.*t13.*1.2852e4+t6.*t14.*1.458e3+signdth1.*t2.*t7.*1.38883248e2)-signdth2.*t6.*3.675e-2,dth1,dth2];
