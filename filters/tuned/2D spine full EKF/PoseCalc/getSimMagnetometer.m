clc

for i=2:1:length(g1.signals.values(:,1))
    th1=camAngles(i,1);
    th2=camAngles(i,2);
    th3=camAngles(i,3);
    th4=camAngles(i,4);


    Rpitch1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
    Rpitch2=[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];
    Rpitch3=[[cos(th3) 0 sin(th3)];[0 1 0];[-sin(th3) 0 cos(th3)]];
    Rpitch4=[[cos(th4) 0 sin(th4)];[0 1 0];[-sin(th4) 0 cos(th4)]];

    R_0_1=(Rpitch1);
    R_0_2=(Rpitch2);
    R_0_3=(Rpitch3);
    R_0_4=(Rpitch4);

    m1temp=R_0_1.'*[magVecX;magVecY;magVecZ];
    m2temp=R_0_2.'*[magVecX;magVecY;magVecZ];
    m3temp=R_0_3.'*[magVecX;magVecY;magVecZ];
    m4temp=R_0_4.'*[magVecX;magVecY;magVecZ];

    m1temp=[m1temp(1);m1temp(2);m1temp(3)]/sqrt(m1temp(1)^2+m1temp(2)^2+m1temp(3)^2);
    m2temp=[m2temp(1);m2temp(2);m2temp(3)]/sqrt(m2temp(1)^2+m2temp(2)^2+m2temp(3)^2);
    m3temp=[m3temp(1);m3temp(2);m3temp(3)]/sqrt(m3temp(1)^2+m3temp(2)^2+m3temp(3)^2);
    m4temp=[m4temp(1);m4temp(2);m4temp(3)]/sqrt(m4temp(1)^2+m4temp(2)^2+m4temp(3)^2);
    
    m1data(i,:)=m1temp;
    m2data(i,:)=m2temp;
    m3data(i,:)=m3temp;
    m4data(i,:)=m4temp;
end

time=0:0.01:t;
ma1Temp.time=time';
ma1Temp.dimensions=3;
ma1Temp.signals.values=m1data;

time=0:0.01:t;
ma2Temp.time=time';
ma2Temp.dimensions=3;
ma2Temp.signals.values=m2data;

time=0:0.01:t;
ma3Temp.time=time';
ma3Temp.dimensions=3;
ma3Temp.signals.values=m3data;

time=0:0.01:t;
ma4Temp.time=time';
ma4Temp.dimensions=3;
ma4Temp.signals.values=m4data;
