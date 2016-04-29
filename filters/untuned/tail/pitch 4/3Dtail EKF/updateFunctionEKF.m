function update = updateFunctionEKF(ddph1,ddph2,ddth1,ddth2,dph1,dph2,dth1,dth2,ph1,ph2,th1,th2)
%UPDATEFUNCTIONEKF
%    UPDATE = UPDATEFUNCTIONEKF(DDPH1,DDPH2,DDTH1,DDTH2,DPH1,DPH2,DTH1,DTH2,PH1,PH2,TH1,TH2)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    21-May-2015 11:54:46

t2 = sin(th1);
t3 = cos(th1);
t4 = cos(ph1);
t5 = sin(ph1);
t6 = dth1.^2;
t7 = ddth1.*t3.*(9.0./2.5e1);
t29 = t2.*t6.*(9.0./2.5e1);
t8 = t7-t29;
t9 = dph1.*t3.*t4.*(9.0./2.5e1);
t34 = dth1.*t2.*t5.*(9.0./2.5e1);
t10 = t9-t34;
t11 = dph1.*t10;
t12 = dph1.*t2.*t5.*(9.0./2.5e1);
t35 = dth1.*t3.*t4.*(9.0./2.5e1);
t13 = t12-t35;
t14 = ddph1.*t3.*t5.*(9.0./2.5e1);
t15 = ddth1.*t2.*t4.*(9.0./2.5e1);
t36 = dth1.*t13;
t16 = t11+t14+t15-t36;
t17 = dph1.*t3.*t5.*(9.0./2.5e1);
t18 = dth1.*t2.*t4.*(9.0./2.5e1);
t19 = t17+t18;
t20 = dph1.*t19;
t21 = dph1.*t2.*t4.*(9.0./2.5e1);
t22 = dth1.*t3.*t5.*(9.0./2.5e1);
t23 = t21+t22;
t24 = dth1.*t23;
t25 = ddth1.*t2.*t5.*(9.0./2.5e1);
t33 = ddph1.*t3.*t4.*(9.0./2.5e1);
t26 = t20+t24+t25-t33;
t27 = sqrt(2.0);
t28 = sin(th2);
t30 = cos(th2);
t31 = sin(ph2);
t32 = cos(ph2);
t37 = dth2.^2;
t38 = ddth2.*t30.*(3.0./1.0e1);
t39 = t7-t29+t38-t28.*t37.*(3.0./1.0e1);
t40 = dph2.*t30.*t32.*(3.0./1.0e1);
t41 = t40-dth2.*t28.*t31.*(3.0./1.0e1);
t42 = dph2.*t41;
t43 = dph2.*t28.*t31.*(3.0./1.0e1);
t44 = t43-dth2.*t30.*t32.*(3.0./1.0e1);
t45 = ddph2.*t30.*t31.*(3.0./1.0e1);
t46 = ddth2.*t28.*t32.*(3.0./1.0e1);
t47 = t11+t14+t15-t36+t42+t45+t46-dth2.*t44;
t48 = dph2.*t30.*t31.*(3.0./1.0e1);
t49 = dth2.*t28.*t32.*(3.0./1.0e1);
t50 = t48+t49;
t51 = dph2.*t50;
t52 = dph2.*t28.*t32.*(3.0./1.0e1);
t53 = dth2.*t30.*t31.*(3.0./1.0e1);
t54 = t52+t53;
t55 = dth2.*t54;
t56 = ddth2.*t28.*t31.*(3.0./1.0e1);
t57 = t20+t24+t25-t33+t51+t55+t56-ddph2.*t30.*t32.*(3.0./1.0e1);
update = [t2.*(-4.9e1./5.0)+t2.*t8-t3.*t4.*t16-t3.*t5.*t26,t3.*(4.9e1./5.0)-t3.*t8-t2.*t4.*t16-t2.*t5.*t26,dth1.*t4,dph1,t3.*t4.*t27.*(1.0./2.0)+t3.*t5.*t27.*(1.0./2.0),t4.*t27.*(1.0./2.0)-t5.*t27.*(1.0./2.0),t2.*t4.*t27.*(1.0./2.0)+t2.*t5.*t27.*(1.0./2.0),t28.*(-4.9e1./5.0)+t28.*t39-t30.*t32.*t47-t30.*t31.*t57,t30.*(4.9e1./5.0)-t30.*t39-t28.*t32.*t47-t28.*t31.*t57,dth2.*t32,dph2,t27.*t30.*t31.*(1.0./2.0)+t27.*t30.*t32.*(1.0./2.0),t27.*t31.*(-1.0./2.0)+t27.*t32.*(1.0./2.0),t27.*t28.*t31.*(1.0./2.0)+t27.*t28.*t32.*(1.0./2.0)];
