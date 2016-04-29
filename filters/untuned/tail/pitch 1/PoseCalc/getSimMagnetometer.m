clc

for i=1:1:length(g1.signals.values(:,1))
    th1=camAngles(i,1);
    th2=camAngles(i,3);
    ph1=camAngles(i,2);
    ph2=camAngles(i,4);


    Rpitch1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
    Ryaw1=[[cos(ph1) -sin(ph1) 0];[sin(ph1) cos(ph1) 0];[0 0 1]];

    Rpitch2=[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];
    Ryaw2=[[cos(ph2) -sin(ph2) 0];[sin(ph2) cos(ph2) 0];[0 0 1]];

    R_0_1=(Ryaw1*Rpitch1);
    R_0_2=(Ryaw2*Rpitch2);


    m1temp=R_0_1.'*[magVecX;magVecY;magVecZ];
    m2temp=R_0_2.'*[magVecX;magVecY;magVecZ];

    m1temp=[m1temp(1);m1temp(2);m1temp(3)]/sqrt(m1temp(1)^2+m1temp(2)^2+m1temp(3)^2);
    m2temp=[m2temp(1);m2temp(2);m2temp(3)]/sqrt(m2temp(1)^2+m2temp(2)^2+m2temp(3)^2);
    
    m1data(i,:)=m1temp;
    m2data(i,:)=m2temp;
end

time=0:0.01:t;
ma1Temp.time=time';
ma1Temp.dimensions=3;
ma1Temp.signals.values=m1data;

time=0:0.01:t;
ma2Temp.time=time';
ma2Temp.dimensions=3;
ma2Temp.signals.values=m2data;
