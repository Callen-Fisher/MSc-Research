clc
clear

%single flick
fileName='data/network.dat';
fd_client1 = fopen('PoseCalc/data/client1.txt');
fd_client2 = fopen('PoseCalc/data/client2.txt');
fd_client3 = fopen('PoseCalc/data/client3.txt');
fd_client4 = fopen('PoseCalc/data/client4.txt'); 
sampleLength=502-5; %must be larger than 105

sampleTimeCam=0.017;%0.0162;
run('PoseCalc/TailFlick3DdataForSim');
run('PoseCalc/calibrateRawDataTailRig');
run('PoseCalc/init')

sim('getDataSim');
run('PoseCalc/getCamAngles');
run('PoseCalc/getSimMagnetometer');
sim('magNoise');
clc

disp('init variables')
disp(' ')
accNoise=0.000213858;
gyroNoise=0.000018992;
magNoise1=0.00000024;%X and Y
magNoise2=0.00000004;%Z
maxA=20;
val1=(maxA*J3)^2*sampleTime^2;
val2=(maxA*J4)^2*sampleTime^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2D EKF
disp('2D EKF')
disp(' ')
addpath('2D spine EKF')
states=[0 0 0 0 0 0 0 0]'; 
covP=diag([1 1 1 1 1 1 1 1]);
I=diag([1 1 1 1 1 1 1 1]);
load('2Dekf')
val3=EKFvalues(1);
val4=EKFvalues(2);
val5=EKFvalues(3);
val6=EKFvalues(4);
clear EKFvalues
R=diag([accNoise,accNoise,gyroNoise,magNoise1,magNoise2,accNoise,accNoise,gyroNoise,magNoise1,magNoise2,accNoise,accNoise,gyroNoise,magNoise1,magNoise2,accNoise,accNoise,gyroNoise,magNoise1,magNoise2]);   
[costAngle2Dekf,costPosition2Dekf,storedStates2Dekf,storedPositions2Dekf]=EKF(val3,val4,val5,val6,states,I,R,a1,a2,a3,a4,sampleTime,g1,g2,g3,g4,covP,camAngles,ma1,ma2,ma3,ma4,l1,l2,l3,l4,camData);
n=length(storedPositions2Dekf(20:end,1))*6;
sumproduct=sum(storedPositions2Dekf(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions2Dekf(20:end,6).*camData.signals.values(20:end,6))+sum(storedPositions2Dekf(20:end,7).*camData.signals.values(20:end,7))+sum(storedPositions2Dekf(20:end,9).*camData.signals.values(20:end,9))+sum(storedPositions2Dekf(20:end,10).*camData.signals.values(20:end,10))+sum(storedPositions2Dekf(20:end,12).*camData.signals.values(20:end,12));
productsum=(sum(storedPositions2Dekf(20:end,4))+sum(storedPositions2Dekf(20:end,6))+sum(storedPositions2Dekf(20:end,7))+sum(storedPositions2Dekf(20:end,9))+sum(storedPositions2Dekf(20:end,10))+sum(storedPositions2Dekf(20:end,12)))*(sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,12)));
sumsquareX=sum(storedPositions2Dekf(20:end,4).^2)+sum(storedPositions2Dekf(20:end,6).^2)+sum(storedPositions2Dekf(20:end,7).^2)+sum(storedPositions2Dekf(20:end,9).^2)+sum(storedPositions2Dekf(20:end,10).^2)+sum(storedPositions2Dekf(20:end,12).^2);
squaresumX=sum(storedPositions2Dekf(20:end,4))+sum(storedPositions2Dekf(20:end,6))+sum(storedPositions2Dekf(20:end,7))+sum(storedPositions2Dekf(20:end,9))+sum(storedPositions2Dekf(20:end,10))+sum(storedPositions2Dekf(20:end,12))^2;
sumsquareY=sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,6).^2)+sum(camData.signals.values(20:end,7).^2)+sum(camData.signals.values(20:end,9).^2)+sum(camData.signals.values(20:end,10).^2)+sum(camData.signals.values(20:end,12).^2);
squaresumY=sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,12))^2;
r2Dekf=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));

disp(['costAngle2Dekf=',num2str(costAngle2Dekf)])
disp(['costPosition2Dekf=',num2str(costPosition2Dekf)])
disp(['r2Dekf=',num2str(r2Dekf)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2D Full EKF
disp(' ')
disp('2D Full EKF')
disp(' ')
addpath('2D spine Full EKF')
states=[0 0 0 0 0 0 0 0 0 0]'; 
covP=diag([1 1 1 1 1 1 1 1 1 1]);
I=diag([1 1 1 1 1 1 1 1 1 1]);
load('2Dfull')
val3=fullEKFvalues(1);
val4=fullEKFvalues(2);
val5=fullEKFvalues(3);
val6=fullEKFvalues(4);
clear fullEKFvalues
R=diag([accNoise,accNoise,gyroNoise,magNoise1,magNoise2,accNoise,accNoise,gyroNoise,magNoise1,magNoise2,accNoise,accNoise,gyroNoise,magNoise1,magNoise2,accNoise,accNoise,gyroNoise,magNoise1,magNoise2]);   
[costAngle2Dfull,costPosition2Dfull,storedStates2Dfull,storedPositions2Dfull]=fullEKF(val1,val2,val3,val4,val5,val6,states,I,R,a1,a2,a3,a4,sampleTime,g1,g2,g3,g4,covP,camAngles,ma1,ma2,ma3,ma4,l1,l2,l3,l4,camData);
n=length(storedPositions2Dfull(20:end,1))*6;
sumproduct=sum(storedPositions2Dfull(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions2Dfull(20:end,6).*camData.signals.values(20:end,6))+sum(storedPositions2Dfull(20:end,7).*camData.signals.values(20:end,7))+sum(storedPositions2Dfull(20:end,9).*camData.signals.values(20:end,9))+sum(storedPositions2Dfull(20:end,10).*camData.signals.values(20:end,10))+sum(storedPositions2Dfull(20:end,12).*camData.signals.values(20:end,12));
productsum=(sum(storedPositions2Dfull(20:end,4))+sum(storedPositions2Dfull(20:end,6))+sum(storedPositions2Dfull(20:end,7))+sum(storedPositions2Dfull(20:end,9))+sum(storedPositions2Dfull(20:end,10))+sum(storedPositions2Dfull(20:end,12)))*(sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,12)));
sumsquareX=sum(storedPositions2Dfull(20:end,4).^2)+sum(storedPositions2Dfull(20:end,6).^2)+sum(storedPositions2Dfull(20:end,7).^2)+sum(storedPositions2Dfull(20:end,9).^2)+sum(storedPositions2Dfull(20:end,10).^2)+sum(storedPositions2Dfull(20:end,12).^2);
squaresumX=sum(storedPositions2Dfull(20:end,4))+sum(storedPositions2Dfull(20:end,6))+sum(storedPositions2Dfull(20:end,7))+sum(storedPositions2Dfull(20:end,9))+sum(storedPositions2Dfull(20:end,10))+sum(storedPositions2Dfull(20:end,12))^2;
sumsquareY=sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,6).^2)+sum(camData.signals.values(20:end,7).^2)+sum(camData.signals.values(20:end,9).^2)+sum(camData.signals.values(20:end,10).^2)+sum(camData.signals.values(20:end,12).^2);
squaresumY=sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,12))^2;
r2Dfull=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));

disp(['costAngle2Dfull=',num2str(costAngle2Dfull)])
disp(['costPosition2Dfull=',num2str(costPosition2Dfull)])
disp(['r2Dfull=',num2str(r2Dfull)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3D EKF
disp(' ')
disp('3D EKF')
disp(' ')
addpath('3D spine EKF')
states=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]'; 
covP=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
I=   diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
load('3Dekf')
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
clear EKFvalues
R=diag([accNoise,accNoise,gyroNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2,...
        accNoise,accNoise,gyroNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2,...
        accNoise,accNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2,...
        accNoise,accNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2]);   
[costPosition3Dekf,storedStates3Dekf,storedPositions3Dekf]=EKF(val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,states,covP,R,I,sampleTime,camData,a1,a2,a3,a4,g1,g2,g3,g4,ma1,ma2,ma3,ma4,l1,l2,l3,l4);
n=length(storedPositions3Dekf(20:end,1))*9;
sumproduct=sum(storedPositions3Dekf(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions3Dekf(20:end,5).*camData.signals.values(20:end,5))+sum(storedPositions3Dekf(20:end,6).*camData.signals.values(20:end,6))+sum(storedPositions3Dekf(20:end,7).*camData.signals.values(20:end,7))+sum(storedPositions3Dekf(20:end,8).*camData.signals.values(20:end,8))+sum(storedPositions3Dekf(20:end,9).*camData.signals.values(20:end,9))+sum(storedPositions3Dekf(20:end,10).*camData.signals.values(20:end,10))+sum(storedPositions3Dekf(20:end,11).*camData.signals.values(20:end,11))+sum(storedPositions3Dekf(20:end,12).*camData.signals.values(20:end,12));
productsum=(sum(storedPositions3Dekf(20:end,4))+sum(storedPositions3Dekf(20:end,5))+sum(storedPositions3Dekf(20:end,6))+sum(storedPositions3Dekf(20:end,7))+sum(storedPositions3Dekf(20:end,8))+sum(storedPositions3Dekf(20:end,9))+sum(storedPositions3Dekf(20:end,10))+sum(storedPositions3Dekf(20:end,11))+sum(storedPositions3Dekf(20:end,12)))*(sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,8))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,11))+sum(camData.signals.values(20:end,12)));
sumsquareX=sum(storedPositions3Dekf(20:end,4).^2)+sum(storedPositions3Dekf(20:end,5).^2)+sum(storedPositions3Dekf(20:end,6).^2)+sum(storedPositions3Dekf(20:end,7).^2)+sum(storedPositions3Dekf(20:end,8).^2)+sum(storedPositions3Dekf(20:end,9).^2)+sum(storedPositions3Dekf(20:end,10).^2)+sum(storedPositions3Dekf(20:end,11).^2)+sum(storedPositions3Dekf(20:end,12).^2);
squaresumX=sum(storedPositions3Dekf(20:end,4))+sum(storedPositions3Dekf(20:end,5))+sum(storedPositions3Dekf(20:end,6))+sum(storedPositions3Dekf(20:end,7))+sum(storedPositions3Dekf(20:end,8))+sum(storedPositions3Dekf(20:end,9))+sum(storedPositions3Dekf(20:end,10))+sum(storedPositions3Dekf(20:end,11))+sum(storedPositions3Dekf(20:end,12))^2;
sumsquareY=sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,5).^2)+sum(camData.signals.values(20:end,6).^2)+sum(camData.signals.values(20:end,7).^2)+sum(camData.signals.values(20:end,8).^2)+sum(camData.signals.values(20:end,9).^2)+sum(camData.signals.values(20:end,10).^2)+sum(camData.signals.values(20:end,11).^2)+sum(camData.signals.values(20:end,12).^2);
squaresumY=sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,8))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,11))+sum(camData.signals.values(20:end,12))^2;
r3Dekf=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));

disp(['costPosition3Dekf=',num2str(costPosition3Dekf)])
disp(['r3Dekf=',num2str(r3Dekf)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3D full EKF
disp(' ')
disp('3D Full EKF')
disp(' ')
addpath('3D spine Full EKF')
covP=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
I=   diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
states=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]'; 
load('3Dfull')
val3=fullEKFvalues(1);
val4=fullEKFvalues(2);
val5=fullEKFvalues(3);
val6=fullEKFvalues(4);
val7=fullEKFvalues(5);
val8=fullEKFvalues(6);
val9=fullEKFvalues(7);
val10=fullEKFvalues(8);
val11=fullEKFvalues(9);
val12=fullEKFvalues(10);
clear fullEKFvalues
R=diag([accNoise,accNoise,gyroNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2,...
        accNoise,accNoise,gyroNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2,...
        accNoise,accNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2,...
        accNoise,accNoise,gyroNoise,gyroNoise,magNoise1,magNoise1,magNoise2]);   
[costPosition3Dfull,storedStates3Dfull,storedPositions3Dfull]=fullEKF(val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,states,covP,R,I,sampleTime,camData,a1,a2,a3,a4,g1,g2,g3,g4,ma1,ma2,ma3,ma4,l1,l2,l3,l4);
n=length(storedPositions3Dfull(20:end,1))*9;
sumproduct=sum(storedPositions3Dfull(20:end,4).*camData.signals.values(20:end,4))+sum(storedPositions3Dfull(20:end,5).*camData.signals.values(20:end,5))+sum(storedPositions3Dfull(20:end,6).*camData.signals.values(20:end,6))+sum(storedPositions3Dfull(20:end,7).*camData.signals.values(20:end,7))+sum(storedPositions3Dfull(20:end,8).*camData.signals.values(20:end,8))+sum(storedPositions3Dfull(20:end,9).*camData.signals.values(20:end,9))+sum(storedPositions3Dfull(20:end,10).*camData.signals.values(20:end,10))+sum(storedPositions3Dfull(20:end,11).*camData.signals.values(20:end,11))+sum(storedPositions3Dfull(20:end,12).*camData.signals.values(20:end,12));
productsum=(sum(storedPositions3Dfull(20:end,4))+sum(storedPositions3Dfull(20:end,5))+sum(storedPositions3Dfull(20:end,6))+sum(storedPositions3Dfull(20:end,7))+sum(storedPositions3Dfull(20:end,8))+sum(storedPositions3Dfull(20:end,9))+sum(storedPositions3Dfull(20:end,10))+sum(storedPositions3Dfull(20:end,11))+sum(storedPositions3Dfull(20:end,12)))*(sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,8))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,11))+sum(camData.signals.values(20:end,12)));
sumsquareX=sum(storedPositions3Dfull(20:end,4).^2)+sum(storedPositions3Dfull(20:end,5).^2)+sum(storedPositions3Dfull(20:end,6).^2)+sum(storedPositions3Dfull(20:end,7).^2)+sum(storedPositions3Dfull(20:end,8).^2)+sum(storedPositions3Dfull(20:end,9).^2)+sum(storedPositions3Dfull(20:end,10).^2)+sum(storedPositions3Dfull(20:end,11).^2)+sum(storedPositions3Dfull(20:end,12).^2);
squaresumX=sum(storedPositions3Dfull(20:end,4))+sum(storedPositions3Dfull(20:end,5))+sum(storedPositions3Dfull(20:end,6))+sum(storedPositions3Dfull(20:end,7))+sum(storedPositions3Dfull(20:end,8))+sum(storedPositions3Dfull(20:end,9))+sum(storedPositions3Dfull(20:end,10))+sum(storedPositions3Dfull(20:end,11))+sum(storedPositions3Dfull(20:end,12))^2;
sumsquareY=sum(camData.signals.values(20:end,4).^2)+sum(camData.signals.values(20:end,5).^2)+sum(camData.signals.values(20:end,6).^2)+sum(camData.signals.values(20:end,7).^2)+sum(camData.signals.values(20:end,8).^2)+sum(camData.signals.values(20:end,9).^2)+sum(camData.signals.values(20:end,10).^2)+sum(camData.signals.values(20:end,11).^2)+sum(camData.signals.values(20:end,12).^2);
squaresumY=sum(camData.signals.values(20:end,4))+sum(camData.signals.values(20:end,5))+sum(camData.signals.values(20:end,6))+sum(camData.signals.values(20:end,7))+sum(camData.signals.values(20:end,8))+sum(camData.signals.values(20:end,9))+sum(camData.signals.values(20:end,10))+sum(camData.signals.values(20:end,11))+sum(camData.signals.values(20:end,12))^2;
r3Dfull=(n*sumproduct-productsum)/(sqrt(n*sumsquareX-squaresumX)*sqrt(n*sumsquareY-squaresumY));


disp(['costPosition3Dfull=',num2str(costPosition3Dfull)])
disp(['r3Dfull=',num2str(r3Dfull)])

run('graphs')