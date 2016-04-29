clc
clear
warning('off')
disp('starting')
fileName='data/network.dat';
fd_client1 = fopen('PoseCalc/data/client1.txt');
fd_client2 = fopen('PoseCalc/data/client2.txt');
fd_client3 = fopen('PoseCalc/data/client3.txt');
fd_client4 = fopen('PoseCalc/data/client4.txt'); 

sampleLength=512-5; %CHANGE
sampleTimeCam=0.0155; %CHANGE 0.145

run('PoseCalc/TailFlick3DdataForSim');
run('PoseCalc/calibrateRawDataTailRig');
run('PoseCalc/init');
sim('getDataSim');
run('PoseCalc/getCamAngles');
run('PoseCalc/getSimMagnetometer');
sim('magNoise');
clc

%%%%%%%%%%%%%%5
disp('initializing noise data')
maxA=20;

accNoise=0.000213858;
gyroNoise=0.000018992;
magNoise_2= 0.0000003435496006870;
magNoise_4= 0.0000002175529393110;
magNoise1=  0.0000004765670605458;
magNoise2=  0.0000000014243194558;
magNoise3=  0.0000006084686349769;
magNoise4=  0.0000000015197751915;
%%%%%%%%%%%%%%%%%%%%%% 2D TRIAD
disp(' ')
disp('2D TRIAD')
disp(' ')
addpath('2D TRIAD')
states=[pi/2 pi/2]'; 
covP=[[1 0];
      [0 1]];
Q = [[gyroNoise  0];
     [0    gyroNoise]];
I=   [[1 0];
      [0 1]];
load('2Dtriad')
val1=TRIADvalues(1);
val2=TRIADvalues(2);
val3=TRIADvalues(3);
val4=TRIADvalues(4);
clear TRIADvalues
[costAngle2Dtriad,costPos2Dtriad,storedStates2Dtriad,storedPositions2Dtriad]=TRIAD(g1,g2,a1,a2,ma1,ma2,magVecX,magVecY,magVecZ,g,Q,I,covP,val1,val2,sampleTime,camAngles,camData,l1,l2,states,val3,val4);
n=length(storedPositions2Dtriad(20:end,1))*6;
sumproduct=sum(storedPositions2Dtriad(20:end,1).*camData.signals.values(20:end,1))+sum(storedPositions2Dtriad(20:end,2).*camData.signals.values(20:end,2))+sum(storedPositions2Dtriad(20:end,3).*camData.signals.values(20:end,3))+sum(storedPositions2Dtriad(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions2Dtriad(20:end,5).*camData.signals.values(20:end,5))+sum(storedPositions2Dtriad(20:end,6).*camData.signals.values(20:end,6));
productsum=(sum(storedPositions2Dtriad(20:end,1))+sum(storedPositions2Dtriad(20:end,2))+sum(storedPositions2Dtriad(20:end,3))+sum(storedPositions2Dtriad(20:end,4))+sum(storedPositions2Dtriad(20:end,5))+sum(storedPositions2Dtriad(20:end,6)))*(sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6)));
sumsquareX=sum(storedPositions2Dtriad(20:end,1).^2)+sum(storedPositions2Dtriad(20:end,2).^2)+sum(storedPositions2Dtriad(20:end,3).^2)+sum(storedPositions2Dtriad(20:end,4).^2)+sum(storedPositions2Dtriad(20:end,5).^2)+sum(storedPositions2Dtriad(20:end,6).^2);
squaresumX=sum(storedPositions2Dtriad(20:end,1))+sum(storedPositions2Dtriad(20:end,2))+sum(storedPositions2Dtriad(20:end,3))+sum(storedPositions2Dtriad(20:end,4))+sum(storedPositions2Dtriad(20:end,5))+sum(storedPositions2Dtriad(20:end,6))^2;
sumsquareY=sum(camData.signals.values(20:end,1).^2)+sum(camData.signals.values(20:end,2).^2)+sum(camData.signals.values(20:end,3).^2)+sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,5).^2)+sum(camData.signals.values(20:end,6).^2);
squaresumY=sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))^2;
r2Dtriad=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));

disp(['costAngle2Dtriad=',num2str(costAngle2Dtriad)])
disp(['costPos2Dtriad=',num2str(costPos2Dtriad)])
disp(['r2Dtriad=',num2str(r2Dtriad)])
%%%%%%%%%%%%%%%%%%%%%% 3D TRIAD
disp(' ')
disp('3D TRIAD')
disp(' ')
addpath('3D TRIAD')
states=[pi/2 0 pi/2 0]'; 
covP=[[1 0 0 0];
      [0 1 0 0];
      [0 0 1 0];
      [0 0 0 1]];
  
I=   [[1 0 0 0];
      [0 1 0 0];
      [0 0 1 0];
      [0 0 0 1]];
load('3Dtriad')
val1=TRIADvalues(1);
val2=TRIADvalues(2);
val3=TRIADvalues(3);
val4=TRIADvalues(4);
val5=TRIADvalues(5);
val6=TRIADvalues(6);
val7=TRIADvalues(7);
val8=TRIADvalues(8);
clear TRIADvalues
Q=diag([gyroNoise  gyroNoise gyroNoise gyroNoise]);
[costAngle3Dtriad,costPos3Dtriad,storedStates3Dtriad,storedPositions3Dtriad]=TRIAD(g1,g2,a1,a2,ma1,ma2,magVecX,magVecY,magVecZ,g,Q,I,covP,val1,val2,val3,val4,sampleTime,camAngles,camData,l1,l2,states,val5,val6,val7,val8);
n=length(storedPositions3Dtriad(20:end,1))*6;
sumproduct=sum(storedPositions3Dtriad(20:end,1).*camData.signals.values(20:end,1))+sum(storedPositions3Dtriad(20:end,2).*camData.signals.values(20:end,2))+sum(storedPositions3Dtriad(20:end,3).*camData.signals.values(20:end,3))+sum(storedPositions3Dtriad(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions3Dtriad(20:end,5).*camData.signals.values(20:end,5))+sum(storedPositions3Dtriad(20:end,6).*camData.signals.values(20:end,6));
productsum=(sum(storedPositions3Dtriad(20:end,1))+sum(storedPositions3Dtriad(20:end,2))+sum(storedPositions3Dtriad(20:end,3))+sum(storedPositions3Dtriad(20:end,4))+sum(storedPositions3Dtriad(20:end,5))+sum(storedPositions3Dtriad(20:end,6)))*(sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6)));
sumsquareX=sum(storedPositions3Dtriad(20:end,1).^2)+sum(storedPositions3Dtriad(20:end,2).^2)+sum(storedPositions3Dtriad(20:end,3).^2)+sum(storedPositions3Dtriad(20:end,4).^2)+sum(storedPositions3Dtriad(20:end,5).^2)+sum(storedPositions3Dtriad(20:end,6).^2);
squaresumX=sum(storedPositions3Dtriad(20:end,1))+sum(storedPositions3Dtriad(20:end,2))+sum(storedPositions3Dtriad(20:end,3))+sum(storedPositions3Dtriad(20:end,4))+sum(storedPositions3Dtriad(20:end,5))+sum(storedPositions3Dtriad(20:end,6))^2;
sumsquareY=sum(camData.signals.values(20:end,1).^2)+sum(camData.signals.values(20:end,2).^2)+sum(camData.signals.values(20:end,3).^2)+sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,5).^2)+sum(camData.signals.values(20:end,6).^2);
squaresumY=sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))^2;
r3Dtriad=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));

disp(['costAngle3Dtriad=',num2str(costAngle3Dtriad)])
disp(['costPos3Dtriad=',num2str(costPos3Dtriad)])
disp(['r3Dtriad=',num2str(r3Dtriad)])
%%%%%%%%%%%%%%%%%%%%%% 2D tail EKF
val1=(maxA*J1)^2*sampleTime^2;
val2=(maxA*J2)^2*sampleTime^2;

disp(' ')
disp('2D tail EKF')
disp(' ')
addpath('2Dtail EKF')
states=[0 0 pi/2 pi/2]'; 
covP=[[1 0 0 0];
      [0 1 0 0];
      [0 0 1 0];
      [0 0 0 1]];
I=   [[1 0 0 0];
      [0 1 0 0];
      [0 0 1 0];
      [0 0 0 1]];
load('2Dekf')
val3=EKFvalues(1);
val4=EKFvalues(2);
clear EKFvalues
R=diag([accNoise,accNoise,gyroNoise,magNoise1,magNoise_2,accNoise,accNoise,gyroNoise,magNoise3,magNoise_4]);  
[costAngle2Dekf,costPos2Dekf,storedStates2Dekf,storedPositions2Dekf]=EKF(val3,val4,states,I,R,a1,a2,sampleTime,g1,g2,covP,camAngles,ma1,ma2,l1,l2,camData);
n=length(storedPositions2Dekf(20:end,1))*4;
sumproduct=sum(storedPositions2Dekf(20:end,1).*camData.signals.values(20:end,1))+sum(storedPositions2Dekf(20:end,3).*camData.signals.values(20:end,3))+sum(storedPositions2Dekf(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions2Dekf(20:end,6).*camData.signals.values(20:end,6));
productsum=(sum(storedPositions2Dekf(20:end,1))+sum(storedPositions2Dekf(20:end,3))+sum(storedPositions2Dekf(20:end,4))+sum(storedPositions2Dekf(20:end,6)))*(sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,6)));
sumsquareX=sum(storedPositions2Dekf(20:end,1).^2)+sum(storedPositions2Dekf(20:end,3).^2)+sum(storedPositions2Dekf(20:end,4).^2)+sum(storedPositions2Dekf(20:end,6).^2);
squaresumX=sum(storedPositions2Dekf(20:end,1))+sum(storedPositions2Dekf(20:end,3))+sum(storedPositions2Dekf(20:end,4))+sum(storedPositions2Dekf(20:end,6))^2;
sumsquareY=sum(camData.signals.values(20:end,1).^2)+sum(camData.signals.values(20:end,3).^2)+sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,6).^2);
squaresumY=sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,6))^2;
r2Dekf=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));


disp(['costAngle2Dekf=',num2str(costAngle2Dekf)])
disp(['costPos2Dekf=',num2str(costPos2Dekf)])
disp(['r2Dekf=',num2str(r2Dekf)])
%%%%%%%%%%%%%%%%%%%%%% 2D tail Full EKF
disp(' ')
disp('2D tail Full EKF')
disp(' ')
addpath('2Dtail full EKF')
states=[0 0 0 0 pi/2 pi/2]'; 
covP=[[1 0 0 0 0 0];
      [0 1 0 0 0 0];
      [0 0 1 0 0 0];
      [0 0 0 1 0 0];
      [0 0 0 0 1 0];
      [0 0 0 0 0 1]];
I=   [[1 0 0 0 0 0];
      [0 1 0 0 0 0];
      [0 0 1 0 0 0];
      [0 0 0 1 0 0];
      [0 0 0 0 1 0];
      [0 0 0 0 0 1]];
load('2Dfull')
val3=fullEKFvalues(1);
val4=fullEKFvalues(2);
clear fullEKFvalues
R=diag([accNoise,accNoise,gyroNoise,magNoise1,magNoise_2,accNoise,accNoise,gyroNoise,magNoise3,magNoise_4]);  
[costAngle2Dfull,costPos2Dfull,storedStates2Dfull,storedPositions2Dfull]=fullEKF(val1,val2,val3,val4,states,I,R,a1,a2,sampleTime,g1,g2,covP,camAngles,ma1,ma2,l1,l2,camData);
n=length(storedPositions2Dfull(20:end,1))*4;
sumproduct=sum(storedPositions2Dfull(20:end,1).*camData.signals.values(20:end,1))+sum(storedPositions2Dfull(20:end,3).*camData.signals.values(20:end,3))+sum(storedPositions2Dfull(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions2Dfull(20:end,6).*camData.signals.values(20:end,6));
productsum=(sum(storedPositions2Dfull(20:end,1))+sum(storedPositions2Dfull(20:end,3))+sum(storedPositions2Dfull(20:end,4))+sum(storedPositions2Dfull(20:end,6)))*(sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,6)));
sumsquareX=sum(storedPositions2Dfull(20:end,1).^2)+sum(storedPositions2Dfull(20:end,3).^2)+sum(storedPositions2Dfull(20:end,4).^2)+sum(storedPositions2Dfull(20:end,6).^2);
squaresumX=sum(storedPositions2Dfull(20:end,1))+sum(storedPositions2Dfull(20:end,3))+sum(storedPositions2Dfull(20:end,4))+sum(storedPositions2Dfull(20:end,6))^2;
sumsquareY=sum(camData.signals.values(20:end,1).^2)+sum(camData.signals.values(20:end,3).^2)+sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,6).^2);
squaresumY=sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,6))^2;
r2Dfull=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));

disp(['costAngle2Dfull=',num2str(costAngle2Dfull)])
disp(['costPos2Dfull=',num2str(costPos2Dfull)])
disp(['r2Dfull=',num2str(r2Dfull)])
%%%%%%%%%%%%%%%%%%%%%% 3D tail EKF
disp(' ')
disp('3D tail EKF')
disp(' ')
addpath('3Dtail EKF')
states=[0 0 0 0 pi/2 0 pi/2 0]'; 
covP=[[1 0 0 0 0 0 0 0];
      [0 1 0 0 0 0 0 0];
      [0 0 1 0 0 0 0 0];
      [0 0 0 1 0 0 0 0];
      [0 0 0 0 1 0 0 0];
      [0 0 0 0 0 1 0 0];
      [0 0 0 0 0 0 1 0];
      [0 0 0 0 0 0 0 1]];
I=   [[1 0 0 0 0 0 0 0];
      [0 1 0 0 0 0 0 0];
      [0 0 1 0 0 0 0 0];
      [0 0 0 1 0 0 0 0];
      [0 0 0 0 1 0 0 0];
      [0 0 0 0 0 1 0 0];
      [0 0 0 0 0 0 1 0];
      [0 0 0 0 0 0 0 1]];
load('3Dekf')
val3=EKFvalues(1);
val4=EKFvalues(2);
val5=EKFvalues(3);
val6=EKFvalues(4);
clear EKFvalues
R=diag([accNoise,accNoise,gyroNoise,gyroNoise,magNoise1,magNoise2,magNoise_2,accNoise,accNoise,gyroNoise,gyroNoise,magNoise3,magNoise4,magNoise_4]);   
[costAngle3Dekf,costPos3Dekf,storedStates3Dekf,storedPositions3Dekf]=EKF(val3,val4,val5,val6,states,I,R,a1,a2,ma1,ma2,sampleTime,g1,g2,covP,camAngles,l1,l2,camData);
n=length(storedPositions3Dekf(20:end,1))*6;
sumproduct=sum(storedPositions3Dekf(20:end,1).*camData.signals.values(20:end,1))+sum(storedPositions3Dekf(20:end,2).*camData.signals.values(20:end,2))+sum(storedPositions3Dekf(20:end,3).*camData.signals.values(20:end,3))+sum(storedPositions3Dekf(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions3Dekf(20:end,5).*camData.signals.values(20:end,5))+sum(storedPositions3Dekf(20:end,6).*camData.signals.values(20:end,6));
productsum=(sum(storedPositions3Dekf(20:end,1))+sum(storedPositions3Dekf(20:end,2))+sum(storedPositions3Dekf(20:end,3))+sum(storedPositions3Dekf(20:end,4))+sum(storedPositions3Dekf(20:end,5))+sum(storedPositions3Dekf(20:end,6)))*(sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6)));
sumsquareX=sum(storedPositions3Dekf(20:end,1).^2)+sum(storedPositions3Dekf(20:end,2).^2)+sum(storedPositions3Dekf(20:end,3).^2)+sum(storedPositions3Dekf(20:end,4).^2)+sum(storedPositions3Dekf(20:end,5).^2)+sum(storedPositions3Dekf(20:end,6).^2);
squaresumX=sum(storedPositions3Dekf(20:end,1))+sum(storedPositions3Dekf(20:end,2))+sum(storedPositions3Dekf(20:end,3))+sum(storedPositions3Dekf(20:end,4))+sum(storedPositions3Dekf(20:end,5))+sum(storedPositions3Dekf(20:end,6))^2;
sumsquareY=sum(camData.signals.values(20:end,1).^2)+sum(camData.signals.values(20:end,2).^2)+sum(camData.signals.values(20:end,3).^2)+sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,5).^2)+sum(camData.signals.values(20:end,6).^2);
squaresumY=sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))^2;
r3Dekf=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));

disp(['costAngle3Dekf=',num2str(costAngle3Dekf)])
disp(['costPos3Dekf=',num2str(costPos3Dekf)])
disp(['r3Dekf=',num2str(r3Dekf)])
%%%%%%%%%%%%%%%%%%%%%% 3D tail Full EKF 
disp(' ')
disp('3D tail Full EKF')
disp(' ')
addpath('3Dtail full EKF')
states=[0 0 0 0 0 0 0 0 pi/2 0 pi/2 0]'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INIT CONDITIONS!!!!!
covP=[[1 0 0 0 0 0 0 0 0 0 0 0];
      [0 1 0 0 0 0 0 0 0 0 0 0];
      [0 0 1 0 0 0 0 0 0 0 0 0];
      [0 0 0 1 0 0 0 0 0 0 0 0];
      [0 0 0 0 1 0 0 0 0 0 0 0];
      [0 0 0 0 0 1 0 0 0 0 0 0];
      [0 0 0 0 0 0 1 0 0 0 0 0];
      [0 0 0 0 0 0 0 1 0 0 0 0];
      [0 0 0 0 0 0 0 0 1 0 0 0];
      [0 0 0 0 0 0 0 0 0 1 0 0];
      [0 0 0 0 0 0 0 0 0 0 1 0];
      [0 0 0 0 0 0 0 0 0 0 0 1]];
I=   [[1 0 0 0 0 0 0 0 0 0 0 0];
      [0 1 0 0 0 0 0 0 0 0 0 0];
      [0 0 1 0 0 0 0 0 0 0 0 0];
      [0 0 0 1 0 0 0 0 0 0 0 0];
      [0 0 0 0 1 0 0 0 0 0 0 0];
      [0 0 0 0 0 1 0 0 0 0 0 0];
      [0 0 0 0 0 0 1 0 0 0 0 0];
      [0 0 0 0 0 0 0 1 0 0 0 0];
      [0 0 0 0 0 0 0 0 1 0 0 0];
      [0 0 0 0 0 0 0 0 0 1 0 0];
      [0 0 0 0 0 0 0 0 0 0 1 0];
      [0 0 0 0 0 0 0 0 0 0 0 1]];
  
load('3Dfull')
val3=fullEKFvalues(1);
val4=fullEKFvalues(2);
val5=fullEKFvalues(3);
val6=fullEKFvalues(4);
clear fullEKFvalues
R=diag([accNoise,accNoise,gyroNoise,gyroNoise,magNoise1,magNoise2,magNoise_2,accNoise,accNoise,gyroNoise,gyroNoise,magNoise3,magNoise4,magNoise_4]); 
[costAngle3Dfull,costPos3Dfull,storedStates3Dfull,storedPositions3Dfull]=fullEKF(val1,val2,val3,val4,val5,val6,states,I,R,a1,a2,ma1,ma2,sampleTime,g1,g2,covP,camAngles,l1,l2,camData);
n=length(storedPositions3Dfull(20:end,1))*6;
sumproduct=sum(storedPositions3Dfull(20:end,1).*camData.signals.values(20:end,1))+sum(storedPositions3Dfull(20:end,2).*camData.signals.values(20:end,2))+sum(storedPositions3Dfull(20:end,3).*camData.signals.values(20:end,3))+sum(storedPositions3Dfull(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions3Dfull(20:end,5).*camData.signals.values(20:end,5))+sum(storedPositions3Dfull(20:end,6).*camData.signals.values(20:end,6));
productsum=(sum(storedPositions3Dfull(20:end,1))+sum(storedPositions3Dfull(20:end,2))+sum(storedPositions3Dfull(20:end,3))+sum(storedPositions3Dfull(20:end,4))+sum(storedPositions3Dfull(20:end,5))+sum(storedPositions3Dfull(20:end,6)))*(sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6)));
sumsquareX=sum(storedPositions3Dfull(20:end,1).^2)+sum(storedPositions3Dfull(20:end,2).^2)+sum(storedPositions3Dfull(20:end,3).^2)+sum(storedPositions3Dfull(20:end,4).^2)+sum(storedPositions3Dfull(20:end,5).^2)+sum(storedPositions3Dfull(20:end,6).^2);
squaresumX=sum(storedPositions3Dfull(20:end,1))+sum(storedPositions3Dfull(20:end,2))+sum(storedPositions3Dfull(20:end,3))+sum(storedPositions3Dfull(20:end,4))+sum(storedPositions3Dfull(20:end,5))+sum(storedPositions3Dfull(20:end,6))^2;
sumsquareY=sum(camData.signals.values(20:end,1).^2)+sum(camData.signals.values(20:end,2).^2)+sum(camData.signals.values(20:end,3).^2)+sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,5).^2)+sum(camData.signals.values(20:end,6).^2);
squaresumY=sum(camData.signals.values(20:end,1))+sum(camData.signals.values(20:end,2))+sum(camData.signals.values(20:end,3))+sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))^2;
r3Dfull=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));

disp(['costAngle3Dfull=',num2str(costAngle3Dfull)])
disp(['costPos3Dfull=',num2str(costPos3Dfull)])
disp(['r3Dfull=',num2str(r3Dfull)])




%run('graphs')






%one beam 3D ekf
disp(' ')
disp('3D tail one beam EKF')
disp(' ')
addpath('oneBeamEKF')
states=[0 0 pi/2 0]'; 
covP=[[1 0 0 0];
      [0 1 0 0];
      [0 0 1 0];
      [0 0 0 1]];
I=   [[1 0 0 0];
      [0 1 0 0];
      [0 0 1 0];
      [0 0 0 1]];

R=diag([accNoise,accNoise,gyroNoise,gyroNoise,magNoise3,magNoise4,magNoise_4]);  

load('oneBeamEKF')
val3=EKFvalues(1);
val4=EKFvalues(2);
clear oneBeamEKF
[costAngle3D,costPos3D,storedStates,storedPositions]=EKF(val3,val4,states,I,R,a1,a2,ma1,ma2,sampleTime,g1,g2,covP,camAngles,l1,l2,camData);
disp(['costPos3D=',num2str(costPos3D)])



%one beam 3D full ekf
disp(' ')
disp('3D tail one beam full EKF')
disp(' ')
addpath('oneBeamFull')
states=[0 0 0 0 pi/2 0]'; 
covP=[[1 0 0 0 0 0];
      [0 1 0 0 0 0];
      [0 0 1 0 0 0];
      [0 0 0 1 0 0];
      [0 0 0 0 1 0];
      [0 0 0 0 0 1]];
  
I=   [[1 0 0 0 0 0];
      [0 1 0 0 0 0];
      [0 0 1 0 0 0];
      [0 0 0 1 0 0];
      [0 0 0 0 1 0];
      [0 0 0 0 0 1]];
R=diag([accNoise,accNoise,gyroNoise,gyroNoise,magNoise3,magNoise4,magNoise_4]);  
load('oneBeamFull')
val3=fullEKFvalues(1);
val4=fullEKFvalues(2);
clear oneBeamFull
[costAngle3DFull,costPos3DFull,storedStates,storedPositions]=fullEKF(val1,val3,val4,states,I,R,a1,a2,ma1,ma2,sampleTime,g1,g2,covP,camAngles,l1,l2,camData);
disp(['costPos3DFull=',num2str(costPos3DFull)])