clc

for i=1:1:length(InertialAngles.signals.values(:,1))
    si1=InertialAngles.signals.values(i,1);
    th1=InertialAngles.signals.values(i,2);
    ph1=InertialAngles.signals.values(i,3);
    
    Rroll1=[[1 0 0];[0 cos(si1) -sin(si1)];[0 sin(si1) cos(si1)]];
    Rpitch1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
    Ryaw1=[[cos(ph1) -sin(ph1) 0];[sin(ph1) cos(ph1) 0];[0 0 1]];

    R=Ryaw1*Rpitch1*Rroll1;

    magCleanData(i,:)=R*[magVecX;magVecY;magVecZ];
    accCleanData(i,:)=R*[0;0;g]+R*InertialAcceleration.signals.values(i,:)';
    gyroCleanData(i,:)=InertialRates.signals.values(i,:)';
end

t=length(magCleanData(:,1))*sampleTime-sampleTime;
time=0:sampleTime:t;
magClean.time=time';
magClean.dimensions=3;
magClean.signals.values=magCleanData;

accClean.time=time';
accClean.dimensions=3;
accClean.signals.values=accCleanData;

gyroClean.time=time';
gyroClean.dimensions=3;
gyroClean.signals.values=gyroCleanData;