clc
clear

%single flick
fileName='Data_ekf8/tail and spine tests/2/2015_05_11_07_23_43.dat';
fd_client1 = fopen('Data_ekf8/tail and spine tests/2/image_data_client1_2015-05-11_07-23.txt');
fd_client2 = fopen('Data_ekf8/tail and spine tests/2/image_data_client2_2015-05-11_07-23.txt');
fd_client3 = fopen('Data_ekf8/tail and spine tests/2/image_data_client3_2015-05-11_07-23.txt');
fd_client4 = fopen('Data_ekf8/tail and spine tests/2/image_data_client4_2015-05-11_07-23.txt'); 
sampleLength=503-5; %must be larger than 105

sampleTimeCam=0.017;%0.0162;
run('TailFlick3DdataForSim');

run('calibrateRawDataTailRig');



J1=10;
J2=10;
J3=0.0114;
J4=0.0114;
sampleTime=0.01;
Cd1=1;
Cd2=1;

rho=1.225;

l1=0.5;
l2=0.3;
l3=0.36;
l4=0.30;

width=0.05;

A1=l3*width;
A2=l4*width;

m1=1.5;
m2=1.5;
m3=0.15;
m4=0.15;

magVecX=1/sqrt(2);
magVecY=1/sqrt(2);
magVecZ=0;

g=9.8;







