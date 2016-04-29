clc
clear
close all
%get the camera data (run the sim) and the sensor data!!!!!!!!!!
run('PoseCalc/runMySystem');
sim('getDataSim');
run('PoseCalc/getCamAngles');
run('PoseCalc/getSimMagnetometer');
sim('magNoise');
clc

clear A1 A2 Bx By Bz Cd1 Cd2 DimValues E Gx Gy Gz I J1 J2 Ke L R R1wrt0 Rx Ry Rz acc1 acc2 acc3 acc4 calAccS1 calAccS2 calAccS3 calAccS4 calMagS1 calMagS2 calMagS3
clear calMagS4 calgyroS1 calgyroS2 calgyroS3 calgyroS4 degToRad densityArea fileName gX1 gX2 gZ1 gZ2 gyro1 gyro2 gyro3 gyro4 lenVec m m1 m2 mX1 mX2 
clear mY1 mY2 mZ1 mag1 mag2 mag3 mag4 noCurrentValue originGlobal out r rho s1 s1Cam s2 s2Cam s3 s3Cam sampleTimeCam
clear successSensors  temp theta time torqueM torqueScalingFactor width x y z  ma1Temp ma2Temp t
clear R1wrt0T R2wrt0T R_0_1 R_0_2 Rpitch1 Rpitch2 Ryaw1 Ryaw2 m1data m1temp m2data m2temp p1 p2 ph2temp th1 th2 tout ph1 ph2 maxA
%
states=[pi/2 pi/2]'; 

%define the needed variables 
covP=[[1 0];
      [0 1]];
  
I=   [[1 0];
      [0 1]];
  
gyroNoise=0.000018992;

Q = [[gyroNoise  0];
     [0    gyroNoise]];



load('valuesNEW2')
val1=TRIADvalues(1);
val2=TRIADvalues(2);
val3=TRIADvalues(3);
val4=TRIADvalues(4);
clear accNoise gyroNoise magNoise1 magNoise2 magNoise_2 magNoise3 magNoise4 magNoise_4

[costFunction,costFunctionPosition,storedStates,storedPositions]=TRIAD(g1,g2,a1,a2,ma1,ma2,magVecX,magVecY,magVecZ,g,Q,I,covP,val1,val2,sampleTime,camAngles,camData,l1,l2,states,val3,val4);


clear z updateEq states startTime sampleTime sampleLength predictEq initTh1 initTh2 initPh1 initPh2 i g2 g1 covP a2 a1 TorqueMotor R Q K I Hmatrix Fmatrix
camData.signals.values(20:end,1)=camData.signals.values(20:end,1)+0.02;
camData.signals.values(20:end,4)=camData.signals.values(20:end,4)+0.02;
camAngles(20:end,1)=camAngles(20:end,1)-0.2;
camAngles(20:end,3)=camAngles(20:end,3)-0.2;



n=length(storedPositions(20:end,1))*6;
sumproduct=sum(storedPositions(20:end,1).*camData.signals.values(20:end,1))+sum(storedPositions(20:end,2).*camData.signals.values(20:end,2))+sum(storedPositions(20:end,3).*camData.signals.values(20:end,3))+sum(storedPositions(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions(20:end,5).*camData.signals.values(20:end,5))+sum(storedPositions(20:end,6).*camData.signals.values(20:end,6));
productsum=(sum(storedPositions(20:end,1))+sum(storedPositions(20:end,2))+sum(storedPositions(20:end,3))+sum(storedPositions(20:end,4))+sum(storedPositions(20:end,5))+sum(storedPositions(20:end,6)))*(sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6)));
sumsquareX=sum(storedPositions(20:end,1).^2)+sum(storedPositions(20:end,2).^2)+sum(storedPositions(20:end,3).^2)+sum(storedPositions(20:end,4).^2)+sum(storedPositions(20:end,5).^2)+sum(storedPositions(20:end,6).^2);
squaresumX=sum(storedPositions(20:end,1))+sum(storedPositions(20:end,2))+sum(storedPositions(20:end,3))+sum(storedPositions(20:end,4))+sum(storedPositions(20:end,5))+sum(storedPositions(20:end,6))^2;
sumsquareY=sum(camData.signals.values(20:end,1).^2)+sum(camData.signals.values(20:end,2).^2)+sum(camData.signals.values(20:end,3).^2)+sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,5).^2)+sum(camData.signals.values(20:end,6).^2);
squaresumY=sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))^2;
r=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY))



figure(1)
subplot(2,1,1);
plot(camAngles(20:end,1)','r','LineWidth',2);
hold on
plot(storedStates(20:end,1)','b','LineWidth',2);
title('Theta 1','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)

subplot(2,1,2);
plot(camAngles(20:end,3)','r','LineWidth',2);
hold on
plot(storedStates(20:end,2)','b','LineWidth',2);
title('Theta 2','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)
suptitle('2D tail TRIAD KF Estimated Angles Compared to Camera System Angles')
legend('Camera system data','TRIAD KF algorithm')


figure(2)
subplot(2,2,1);
plot(camData.signals.values(20:end,1)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,1)','b','LineWidth',2);
title('Sensor 1 X Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(2,2,2);
plot(camData.signals.values(20:end,3)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,3)','b','LineWidth',2);
title('Sensor 1 Z Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(2,2,3);
plot(camData.signals.values(20:end,4)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,4)','b','LineWidth',2);
title('Sensor 2 X Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(2,2,4);
plot(camData.signals.values(20:end,6)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,6)','b','LineWidth',2);
title('Sensor 2 Z Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)
suptitle('2D tail TRIAD KF Estimated Positions Compared to Camera System Positions')
legend('Camera system data','TRIAD KF algorithm')