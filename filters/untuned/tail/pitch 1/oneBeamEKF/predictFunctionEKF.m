function predict = predictFunctionEKF(dph,dth,ph,th)
%PREDICTFUNCTIONEKF
%    PREDICT = PREDICTFUNCTIONEKF(DPH,DTH,PH,TH)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    02-Jul-2015 12:08:10

t2 = ph.*2.0;
t3 = sin(th);
t4 = sin(t2);
t5 = th.*2.0;
t6 = sin(t5);
predict = [(dph.^2.*t6.*(-5.445e2)+dph.*dth.*t4.*1.52e3)./(cos(t2).*7.6e2+1.849e3),-(dth.*(dph.*t6.*1.089e3-dth.*t4.*7.6e2))./(t3.^2.*1.089e3-2.609e3),dth,dph];