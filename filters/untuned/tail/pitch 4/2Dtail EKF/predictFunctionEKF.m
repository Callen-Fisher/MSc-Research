function predict = predictFunctionEKF(dth1,dth2,th1,th2)
%PREDICTFUNCTIONEKF
%    PREDICT = PREDICTFUNCTIONEKF(DTH1,DTH2,TH1,TH2)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    23-May-2015 10:56:00

t3 = th1-th2;
t2 = sin(t3);
t4 = t2.^2;
t5 = t4.*2.916e3;
t6 = t5+2.0527e4;
t7 = 1.0./t6;
t8 = dth1.^2;
t9 = dth2.^2;
t10 = th1.*2.0;
t11 = t10-th2.*2.0;
t12 = sin(t11);
predict = [t7.*(t2.*t9.*1.97e2+t8.*t12.*5.4e1).*-2.7e1,t7.*(t2.*t8.*1.19e2+t9.*t12.*(2.7e1./2.0)).*1.08e2,dth1,dth2];