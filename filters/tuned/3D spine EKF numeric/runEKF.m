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

maxA=20;
val1=(maxA*J3)^2*sampleTime^2;
val2=(maxA*J4)^2*sampleTime^2;

clear A1 A2 Bx By Bz Cd1 Cd2 DimValues E Gx Gy Gz I J1 J2 J3 J4 Ke L R R1wrt0 Rx Ry Rz acc1 acc2 acc3 acc4 calAccS1 calAccS2 calAccS3 calAccS4 calMagS1 calMagS2 calMagS3
clear calMagS4 calgyroS1 calgyroS2 calgyroS3 calgyroS4 degToRad densityArea fileName gX1 gX2 gZ1 gZ2 gyro1 gyro2 gyro3 gyro4 lenVec m m1 m2 m3 m4 mX1 mX2 
clear mY1 mY2 mZ1 mag1 mag2 mag3 mag4 magVecX magVecY magVecZ noCurrentValue originGlobal out r rho s1 s1Cam s2 s2Cam s3 s3Cam sampleTimeCam
clear successSensors  temp theta time torqueM torqueScalingFactor width x y z g ma1Temp ma2Temp t
clear R1wrt0T R2wrt0T R_0_1 R_0_2 R_0_3 R_0_4 Rpitch1 Rpitch2 Rpitch3 Rpitch4 Ryaw1 Ryaw2 Ryaw3 Ryaw4 m1data m1temp m2data m2temp m3data m3temp m4data m4temp p1 p2 ph2temp th1 th2 tout ph1 ph2 maxA
clear i ma3Temp ma4Temp th3 th4 

states=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]'; 

%define the needed variables 
covP=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
  
I=   diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
  

accNoise=0.000213858;
gyroNoise=0.000018992;

magNoise1=0.00000024;%X and Y
magNoise2=0.00000004;%Z

load('values')
val3=EKFvalues(1);
val4=EKFvalues(2);
val5=EKFvalues(3);
val6=EKFvalues(4);
val7=EKFvalues(5);
val8=EKFvalues(6);
val9=EKFvalues(7);
val10=EKFvalues(8);
val11=EKFvalues(9);
val12=EKFvalues(10);


R=diag([accNoise,accNoise,gyroNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2,...
        accNoise,accNoise,gyroNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2,...
        accNoise,accNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2,...
        accNoise,accNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2]);   
clear accNoise gyroNoise magNoise1 magNoise2 magNoise_2 magNoise3 magNoise4 magNoise_4

[costFunction,storedStates,storedPositions]=EKF(val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,states,covP,R,I,sampleTime,camData,a1,a2,a3,a4,g1,g2,g3,g4,ma1,ma2,ma3,ma4,l1,l2,l3,l4);

clear z updateEq states startTime sampleTime sampleLength predictEq initTh1 initTh2 initPh1 initPh2 i g2 g1 covP a2 a1 TorqueMotor R Q K I Hmatrix Fmatrix



n=length(storedPositions(20:end,1))*9;
sumproduct=sum(storedPositions(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions(20:end,5).*camData.signals.values(20:end,5))+sum(storedPositions(20:end,6).*camData.signals.values(20:end,6))+sum(storedPositions(20:end,7).*camData.signals.values(20:end,7))+sum(storedPositions(20:end,8).*camData.signals.values(20:end,8))+sum(storedPositions(20:end,9).*camData.signals.values(20:end,9))+sum(storedPositions(20:end,10).*camData.signals.values(20:end,10))+sum(storedPositions(20:end,11).*camData.signals.values(20:end,11))+sum(storedPositions(20:end,12).*camData.signals.values(20:end,12));
productsum=(sum(storedPositions(20:end,4))+sum(storedPositions(20:end,5))+sum(storedPositions(20:end,6))+sum(storedPositions(20:end,7))+sum(storedPositions(20:end,8))+sum(storedPositions(20:end,9))+sum(storedPositions(20:end,10))+sum(storedPositions(20:end,11))+sum(storedPositions(20:end,12)))*(sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,8))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,11))+sum(camData.signals.values(20:end,12)));
sumsquareX=sum(storedPositions(20:end,4).^2)+sum(storedPositions(20:end,5).^2)+sum(storedPositions(20:end,6).^2)+sum(storedPositions(20:end,7).^2)+sum(storedPositions(20:end,8).^2)+sum(storedPositions(20:end,9).^2)+sum(storedPositions(20:end,10).^2)+sum(storedPositions(20:end,11).^2)+sum(storedPositions(20:end,12).^2);
squaresumX=sum(storedPositions(20:end,4))+sum(storedPositions(20:end,5))+sum(storedPositions(20:end,6))+sum(storedPositions(20:end,7))+sum(storedPositions(20:end,8))+sum(storedPositions(20:end,9))+sum(storedPositions(20:end,10))+sum(storedPositions(20:end,11))+sum(storedPositions(20:end,12))^2;
sumsquareY=sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,5).^2)+sum(camData.signals.values(20:end,6).^2)+sum(camData.signals.values(20:end,7).^2)+sum(camData.signals.values(20:end,8).^2)+sum(camData.signals.values(20:end,9).^2)+sum(camData.signals.values(20:end,10).^2)+sum(camData.signals.values(20:end,11).^2)+sum(camData.signals.values(20:end,12).^2);
squaresumY=sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,8))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,11))+sum(camData.signals.values(20:end,12))^2;
r=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY))


close all
figure(1)
clf;
subplot(3,3,1);
plot(camData.signals.values(20:end,4)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,4)','b','LineWidth',2);
title('Sensor 2 X Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(3,3,2);
plot(camData.signals.values(20:end,5)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,5)','b','LineWidth',2);
title('Sensor 2 Y Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(3,3,3);
plot(camData.signals.values(20:end,6)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,6)','b','LineWidth',2);
title('Sensor 2 Z Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(3,3,4);
plot(camData.signals.values(20:end,7)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,7)','b','LineWidth',2);
title('Sensor 3 X Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(3,3,5);
plot(camData.signals.values(20:end,8)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,8)','b','LineWidth',2);
title('Sensor 3 Y Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(3,3,6);
plot(camData.signals.values(20:end,9)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,9)','b','LineWidth',2);
title('Sensor 3 Z Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(3,3,7);
plot(camData.signals.values(20:end,10)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,10)','b','LineWidth',2);
title('Sensor 4 X Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(3,3,8);
plot(camData.signals.values(20:end,11)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,11)','b','LineWidth',2);
title('Sensor 4 Y Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

subplot(3,3,9);
plot(camData.signals.values(20:end,12)','r','LineWidth',2);
hold on
plot(storedPositions(20:end,12)','b','LineWidth',2);
title('Sensor 4 Z Position','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Position (m)','FontSize',10)

suptitle('3D Spine and Tail Multibody Estimated Positions Compared to Camera System Positions')
legend('Camera system data','Multibody algorithm')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








figure(2)
clf;
subplot(4,3,1);
plot(storedStates(20:end,11)','b','LineWidth',2);
title('Phi 1','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)
subplot(4,3,2);
plot(storedStates(20:end,12)','b','LineWidth',2);
title('Theta 1','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)
subplot(4,3,3);
plot(storedStates(20:end,13)','b','LineWidth',2);
title('Psi 1','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)

subplot(4,3,4);
plot(storedStates(20:end,14)','b','LineWidth',2);
title('Phi 2','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)
subplot(4,3,5);
plot(storedStates(20:end,15)','b','LineWidth',2);
title('Theta 2','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)
subplot(4,3,6);
plot(storedStates(20:end,16)','b','LineWidth',2);
title('Psi 2','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)


subplot(4,3,8);
plot(storedStates(20:end,17)','b','LineWidth',2);
title('Theta 3','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)
subplot(4,3,9);
plot(storedStates(20:end,18)','b','LineWidth',2);
title('Psi 3','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)


subplot(4,3,11);
plot(storedStates(20:end,19)','b','LineWidth',2);
title('Theta 4','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)
subplot(4,3,12);
plot(storedStates(20:end,20)','b','LineWidth',2);
title('Psi 4','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angle (radians)','FontSize',10)
suptitle('3D Spine and Tail Multibody Estimated Angles')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
clf;
subplot(4,3,1);
plot(storedStates(20:end,1)','b','LineWidth',2);
title('Phi 1 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)
subplot(4,3,2);
plot(storedStates(20:end,2)','b','LineWidth',2);
title('Theta 1 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)
subplot(4,3,3);
plot(storedStates(20:end,3)','b','LineWidth',2);
title('Psi 1 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)

subplot(4,3,4);
plot(storedStates(20:end,4)','b','LineWidth',2);
title('Phi 2 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)
subplot(4,3,5);
plot(storedStates(20:end,5)','b','LineWidth',2);
title('Theta 2 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)
subplot(4,3,6);
plot(storedStates(20:end,6)','b','LineWidth',2);
title('Psi 2 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)

subplot(4,3,8);
plot(storedStates(20:end,7)','b','LineWidth',2);
title('Theta 3 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)
subplot(4,3,9);
plot(storedStates(20:end,8)','b','LineWidth',2);
title('Psi 3 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)

subplot(4,3,11);
plot(storedStates(20:end,9)','b','LineWidth',2);
title('Theta 4 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)
subplot(4,3,12);
plot(storedStates(20:end,10)','b','LineWidth',2);
title('Psi 4 angular rate estimate','FontSize',12)
xlabel('Samples','FontSize',10)
ylabel('Angular Rate (rad/s)','FontSize',10)

suptitle('3D Spine and Tail Multibody Estimated Angular Rate')