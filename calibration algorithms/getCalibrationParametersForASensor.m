clc 
clear

%the file names:

accXdown='accXdown.dat'; %positiove 
accXup='accXup.dat';
accYdown='accYdown.dat';
accYup='accYup.dat';
accZdown='accZdown.dat';
accZup='accZup.dat';

gyroData='gyro.dat';
magData='mag.dat';

%calibrate the accelerometer
calibrateAccFunction(accZdown,accZup,accYdown,accYup,accXdown,accXup);
%calibrate the magnetometer
calibrateMagFunction(magData);
%calibrate the gyro
calibrateGyroFunction(gyroData);