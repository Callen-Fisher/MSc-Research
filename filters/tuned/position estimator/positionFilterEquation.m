clc;
clear;
disp('starting');
syms si th ph dsi dth dph ax ay az vx vy vz px py pz g magVecX magVecY magVecZ
syms accX accY accZ gyroX gyroY gyroZ

Rr=[[1 0 0];[0 cos(si) -sin(si)];[0 sin(si) cos(si)]];
Rp=[[cos(th) 0 sin(th)];[0 1 0];[-sin(th) 0 cos(th)]];
Ry=[[cos(ph) -sin(ph) 0];[sin(ph) cos(ph) 0];[0 0 1]];

R=Ry*Rp*Rr;
R_gyro=[[cos(ph)*cos(th) -sin(ph) 0];[sin(ph)*cos(th) cos(ph) 0];[-sin(th) 0 1]]; %%%%INERTIA TO BODY


disp('the Kalman Filter code');

disp('the F matrix');
states=[si th ph vx vy vz px py pz];

temp=R*[accX,accY,accZ]'-[0;0;g];
ax=temp(1);
ay=temp(2);
az=temp(3);
     
R_gyro_B_I=[[1 0 sin(th)];[0 cos(si) -sin(si)*cos(th)];[0 sin(si) cos(si)*cos(th)]];
temp=R_gyro_B_I*[gyroX,gyroY,gyroZ]';
gx=temp(1);
gy=temp(2);
gz=temp(3);

Fequations=[gx gy gz ax ay az vx vy vz];
F = jacobian(Fequations,states);

disp('the measurement equations')
acc=R.'*[0;0;g]+R.'*[ax;ay;az];
ax=acc(1);
ay=acc(2);
az=acc(3);
gyro=R_gyro*[dsi;dth;dph];
gyroX=gyro(1);
gyroY=gyro(2);
gyroZ=gyro(3);
mag=R.'*[magVecX;magVecY;magVecZ];
magX=mag(1);
magY=mag(2);
magZ=mag(3);
gps=[px py pz];
gpsX=gps(1);
gpsY=gps(2);
gpsZ=gps(3);

disp('the H matrix');
Hequations1=[ax ay az gyroX gyroY gyroZ magX magY magZ];
Hequations2=[ax ay az gyroX gyroY gyroZ magX magY magZ gpsX gpsY gpsZ];

H1 = jacobian(Hequations1,states);
H2 = jacobian(Hequations2,states);

predict=Fequations;
update1=Hequations1;
update2=Hequations2;

disp('matlab function F')
matlabFunction(F,'file','Ffunc');
disp('matlab function H')
matlabFunction(H1,'file','Hfunc1');
matlabFunction(H2,'file','Hfunc2');
disp('matlab function predict')
matlabFunction(predict,'file','predictFunc');
disp('matlab function update')
matlabFunction(update1,'file','updateFunc1');
matlabFunction(update2,'file','updateFunc2');

clear all
clc