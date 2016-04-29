function [ U,c ] = calibrateMagFunction( fileName )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%the file format:
%GyroX1,GyroY1,GyroZ1,AccX1,AccY1,AccZ1,MagX1,MagY1,MagZ1,Temp1,rollCommand,pitchCommand,yawCommand,rollEncoder,pitchEncoder,yawEncoder,yawVel

%the file names:
%Mag.dat

[ gyroS1, magS1, accS1, tempS1, tempADCS1, yawEncoder,yawVelocity] = calibrationReadDataFunction( fileName );
save('rawMag','magS1')
[U,c] = MgnCalibration(magS1)
disp('U     :  shape ellipsoid parameter, (3x3) upper triangular matrix');
disp('c      : ellipsoid center, (3x1) vector');
disp('');
disp('Ellipsoid equation : transpose(v-c)*(transpose(U)*U)(v-c) = 1 ');
disp('with v a rough triaxes magnetometer  measurement');
disp('');
disp('calibrated measurement w = U*(v-c)');

end

