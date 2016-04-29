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
clear mY1 mY2 mZ1 mag1 mag2 mag3 mag4 magVecX magVecY magVecZ noCurrentValue originGlobal out r rho s1 s1Cam s2 s2Cam s3 s3Cam sampleTimeCam
clear successSensors  temp theta time torqueM torqueScalingFactor width x y z g ma1Temp ma2Temp t
clear R1wrt0T R2wrt0T R_0_1 R_0_2 Rpitch1 Rpitch2 Ryaw1 Ryaw2 m1data m1temp m2data m2temp p1 p2 ph2temp th1 th2 tout ph1 ph2 maxA

states=[0 0 pi/2 0]'; 

%define the needed variables 
covP=[[1 0 0 0];
      [0 1 0 0];
      [0 0 1 0];
      [0 0 0 1]];
  
I=   [[1 0 0 0];
      [0 1 0 0];
      [0 0 1 0];
      [0 0 0 1]];
  

accNoise=0.000213858;
gyroNoise=0.000018992;

magNoise_2= 0.0000003435496006870;%added 4 zeros
magNoise_4= 0.0000002175529393110;
magNoise1=  0.0000004765670605458;
magNoise2=  0.0000000014243194558;
magNoise3=  0.0000006084686349769;
magNoise4=  0.0000000015197751915;

R=diag([accNoise,accNoise,gyroNoise,gyroNoise,magNoise3,magNoise4,magNoise_4]);  

load('values')
val3=EKFvalues(1);
val4=EKFvalues(2);


 
clear accNoise gyroNoise magNoise1 magNoise2 magNoise_2 magNoise3 magNoise4 magNoise_4

[costFunction,costFunctionPosition,storedStates,storedPositions]=EKF(val3,val4,states,I,R,a1,a2,ma1,ma2,sampleTime,g1,g2,covP,camAngles,l1,l2,camData);


clear z updateEq states startTime sampleTime sampleLength predictEq initTh1 initTh2 initPh1 initPh2 i g2 g1 covP a2 a1 TorqueMotor R Q K I Hmatrix Fmatrix


subplot(1,3,1);
plot(camData.signals.values(20:end,4)','r');
hold on
plot(storedPositions(20:end,1)','b');
title('Sensor 1 X Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(1,3,2);
plot(camData.signals.values(20:end,5)','r');
hold on
plot(storedPositions(20:end,2)','b');
title('Sensor 1 Y Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(1,3,3);
plot(camData.signals.values(20:end,6)','r');
hold on
plot(storedPositions(20:end,3)','b');
title('Sensor 1 Z Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)
suptitle('3D tail Multibody Estimated Positions Compared to Camera System Positions')
legend('Camera system data','Multibody algorithm')