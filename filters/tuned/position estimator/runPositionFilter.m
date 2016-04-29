clc;
clear;
magVecX=sqrt(2)/2;
magVecY=sqrt(2)/2;
magVecZ=0;
g=9.81;
gpsSampleTime=0.1;
sampleTime=0.01;
sim('simulatedData');
run('simSensors');
sim('addSensorNoise');

accNoise=0.001;
magNoise=0.001;
gyroNoise=0.001;
gpsNoise=0.001;
states=[0,0,0,20,0,0,0,0,0]';
covP=diag([1,1,1,1,1,1,1,1,1]);
I=diag([1,1,1,1,1,1,1,1,1]);
R1=diag([accNoise,accNoise,accNoise,gyroNoise,gyroNoise,gyroNoise,magNoise,magNoise,magNoise]);
R2=diag([accNoise,accNoise,accNoise,gyroNoise,gyroNoise,gyroNoise,magNoise,magNoise,magNoise,gpsNoise,gpsNoise,gpsNoise]);

load('values')
val1=POSITIONvalues(1);
val2=POSITIONvalues(2);
val3=POSITIONvalues(3);
val4=POSITIONvalues(4);
val5=POSITIONvalues(5);
val6=POSITIONvalues(6);

[costFunction,storedStates,storedPositions]=positionFilter(acc,mag,gyro,g,magVecX,magVecY,magVecZ,R1,R2,I,covP,val1,val2,val3,val4,val5,val6,truePosition,states,sampleTime,GPS_position);

figure(1)
subplot(1,3,1);
plot(truePosition.signals.values(:,1)','r','LineWidth',2);
hold on
plot(storedPositions(:,1)','b','LineWidth',2);
title('X Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(1,3,2);
plot(truePosition.signals.values(:,2)','r','LineWidth',2);
hold on
plot(storedPositions(:,2)','b','LineWidth',2);
title('Y Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(1,3,3);
plot(truePosition.signals.values(:,3)','r','LineWidth',2);
hold on
plot(storedPositions(:,3)','b','LineWidth',2);
title('Z Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)
suptitle('Position Estimate with No GPS Update')
legend('True Position','Estimated Position')