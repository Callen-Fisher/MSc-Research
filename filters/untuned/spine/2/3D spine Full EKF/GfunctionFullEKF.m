function G = GfunctionFullEKF(th1,th2,th3,th4)
%GFUNCTIONFULLEKF
%    G = GFUNCTIONFULLEKF(TH1,TH2,TH3,TH4)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    13-Jun-2015 12:56:50

G = [0.0;cos(th1).*(-1.47e2./1.0e2);0.0;0.0;cos(th2).*(-4.41e2./5.0e2);0.0;cos(th3).*(-7.938e-1);0.0;cos(th4).*(-2.205e-1);0.0];
