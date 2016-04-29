function update1 = updateFunc1(accX,accY,accZ,dph,dsi,dth,g,magVecX,magVecY,magVecZ,ph,si,th)
%UPDATEFUNC1
%    UPDATE1 = UPDATEFUNC1(ACCX,ACCY,ACCZ,DPH,DSI,DTH,G,MAGVECX,MAGVECY,MAGVECZ,PH,SI,TH)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    03-Jun-2015 17:23:21

t2 = sin(th);
t3 = cos(th);
t4 = conj(accY);
t5 = cos(si);
t6 = cos(ph);
t7 = sin(si);
t8 = conj(accZ);
t9 = sin(ph);
t10 = conj(accX);
t11 = t5.*t9;
t13 = t2.*t6.*t7;
t12 = t11-t13;
t14 = t7.*t9;
t15 = t2.*t5.*t6;
t16 = t14+t15;
t17 = t8.*t16;
t18 = t3.*t6.*t10;
t30 = t4.*t12;
t19 = t17+t18-t30;
t20 = t5.*t6;
t21 = t2.*t7.*t9;
t22 = t20+t21;
t23 = t4.*t22;
t24 = t6.*t7;
t31 = t2.*t5.*t9;
t25 = t24-t31;
t26 = t3.*t9.*t10;
t32 = t8.*t25;
t27 = t23+t26-t32;
t28 = t2.*t10;
t33 = t3.*t5.*t8;
t34 = t3.*t4.*t7;
t29 = g+t28-t33-t34;
update1 = [-g.*t2+t2.*t29+t3.*t6.*t19+t3.*t9.*t27,-t12.*t19+t22.*t27+g.*t3.*t7-t3.*t7.*t29,t16.*t19-t25.*t27+g.*t3.*t5-t3.*t5.*t29,-dth.*t9+dsi.*t3.*t6,dth.*t6+dsi.*t3.*t9,dph-dsi.*t2,-magVecZ.*t2+magVecX.*t3.*t6+magVecY.*t3.*t9,-magVecX.*t12+magVecY.*t22+magVecZ.*t3.*t7,magVecX.*t16-magVecY.*t25+magVecZ.*t3.*t5];
