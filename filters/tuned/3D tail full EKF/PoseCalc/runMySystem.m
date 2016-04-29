clc
clear

%single flick
fileName='Data_ekf8/singleTest/2015_02_25_13_21_56.dat';
fd_client1 = fopen('Data_ekf8/singleTest/image_data_client1_2015-02-25_13-18.txt');
fd_client2 = fopen('Data_ekf8/singleTest/image_data_client2_2015-02-25_13-18.txt');
fd_client3 = fopen('Data_ekf8/singleTest/image_data_client3_2015-02-25_13-18.txt');
fd_client4 = fopen('Data_ekf8/singleTest/image_data_client4_2015-02-25_13-18.txt'); 
sampleLength=399-5; %must be larger than 105

sampleTimeCam=0.0145;
run('TailFlick3DdataForSim');

run('calibrateRawDataTailRig');

mY1=mean(mag2.signals.values(1:40,2));
mX1=mean(mag2.signals.values(1:40,1));
mY2=mean(mag3.signals.values(1:40,2));
mX2=mean(mag3.signals.values(1:40,1));
mZ1=mean(mag2.signals.values(1:40,3));

gX1=mean(acc2.signals.values(1:40,1));
gZ1=mean(acc2.signals.values(1:40,3));
gX2=mean(acc3.signals.values(1:40,1));
gZ2=mean(acc3.signals.values(1:40,3));

initPh1=atan2(-mY1,mX1);
initPh2=atan2(-mY2,mX2)-initPh1;
initTh1=atan2(-gX1,gZ1);
initTh2=atan2(-gX2,gZ2)-initTh1;

J1=0.0114;
J2=0.0114;

sampleTime=0.01;
Cd1=1;
Cd2=1;

rho=1.225;

l1=0.36;
l2=0.30;

width=0.05;

A1=l1*width;
A2=l2*width;

m1=0.15;
m2=0.15;

initPh1=0*pi/180;
R1wrt0=[[cos(initPh1)*cos(initTh1) -sin(initPh1)    cos(initPh1)*sin(initTh1)];...
        [ sin(initPh1)*cos(initTh1)             cos(initPh1)                 sin(initPh1)*sin(initTh1) ];...
        [ -sin(initTh1)  0  cos(initTh1)]]';
    
temp=R1wrt0'*[mX1,mY1,mZ1]';
lenVec=sqrt(mZ1^2+mY1^2+mX1^2);
magVecX=temp(1)/lenVec;
magVecY=temp(2)/lenVec;
magVecZ=temp(3)/lenVec;

magVecY=-magVecY;
magVecZ=2*magVecZ;

lenVec=sqrt(magVecX^2+magVecY^2+magVecZ^2);
magVecX=magVecX/lenVec;
magVecY=0.5*magVecY/lenVec;
magVecZ=magVecZ/lenVec;

magVecX=1/sqrt(2);
magVecY=1/sqrt(2);
magVecZ=0;

g=9.8;




%%%%NEW SIM
L=l1+l2;
m=0.300;
r=0.05;
densityArea=m/(pi*r^2*L);
E=0.1;
I=1/4*m*r^2+1/3*m*L^2;
Ke=0.1;






