function [ X ] = calibrateAccFunction( ZdownFileName,ZupFileName,YdownFileName,YupFileName,XdownFileName,XupFileName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%the file format:
%GyroX1,GyroY1,GyroZ1,AccX1,AccY1,AccZ1,MagX1,MagY1,MagZ1,Temp1,Temp2,rollCommand,pitchCommand,yawCommand,rollEncoder,pitchEncoder,yawEncoder,yawVel

%the file names:
%Zdown, Zup, Ydown, Yup, Xdown, Xup

%calibrated=[raw,1]*A
%therefore Y=W.X    X is the 12 calibration values
%Y= known and normalized gravity vectors
%w= sensor raw data
%X= the 12 calibration values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POSITION 1
%Zb down: [Ax1 Ay1 Az1]=[0 0 1] 
[gyroS1, magS1, accZDown, tempS1, tempADCS1, yawEncoder,yawVelocity]=calibrationReadDataFunction(ZdownFileName);
disp('Zdown file complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POSITION 2
%Zb up: [Ax1 Ay1 Az1]=[0 0 -1] 
[gyroS1, magS1, accZUp, tempS1, tempADCS1, yawEncoder,yawVelocity]=calibrationReadDataFunction(ZupFileName);
disp('Zup file complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POSITION 3
%Yb down: [Ax1 Ay1 Az1]=[0 1 0] 
[gyroS1, magS1, accYDown, tempS1, tempADCS1, yawEncoder,yawVelocity]=calibrationReadDataFunction(YdownFileName);
disp('Ydown file complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POSITION 4
%Yb up: [Ax1 Ay1 Az1]=[0 -1 0]
[gyroS1, magS1, accYUp, tempS1, tempADCS1, yawEncoder,yawVelocity]=calibrationReadDataFunction(YupFileName);
disp('Yup file complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POSITION 5
%Xb down: [Ax1 Ay1 Az1]=[1 0 0] 
[gyroS1, magS1, accXDown, tempS1, tempADCS1, yawEncoder,yawVelocity]=calibrationReadDataFunction(XdownFileName);
disp('Xdown file complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POSITION 6
%Xb down: [Ax1 Ay1 Az1]=[-1 0 0] 
[gyroS1, magS1, accXUp, tempS1, tempADCS1, yawEncoder,yawVelocity]=calibrationReadDataFunction(XupFileName);
disp('Xup file complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%combine them 
Y1=[];
w1=[];
for i=1:1:length(accZDown)
    Y1=[Y1;[0 0 1]];
    w1=[w1;[accZDown(i,:) 1]];
end
Y2=[];
w2=[];
for i=1:1:length(accZUp)
    Y2=[Y2;[0 0 -1]];
    w2=[w2;[accZUp(i,:) 1]];
end

Y3=[];
w3=[];
for i=1:1:length(accYDown)
    Y3=[Y3;[0 1 0]];
    w3=[w3;[accYDown(i,:) 1]];
end
Y4=[];
w4=[];
for i=1:1:length(accYUp)
    Y4=[Y4;[0 -1 0]];
    w4=[w4;[accYUp(i,:) 1]];
end

Y5=[];
w5=[];
for i=1:1:length(accXDown)
    Y5=[Y5;[1 0 0]];
    w5=[w5;[accXDown(i,:) 1]];
end
Y6=[];
w6=[];
for i=1:1:length(accXUp)
    Y6=[Y6;[-1 0 0]];
    w6=[w6;[accXUp(i,:) 1]];
end
Y=[[Y1]
   [Y2]
   [Y3]
   [Y4]
   [Y5]
   [Y6]];
W=[[w1]
   [w2]
   [w3]
   [w4]
   [w5]
   [w6]];

%calculate the 3 bias values:
biasZ=mean(accZUp(:,3))+1;
biasY=mean(accYUp(:,2))+1;
biasX=mean(accXUp(:,1))+1;

%solve for the 12 unknowns 
X=[[];[]];
X=inv(transpose(W)*W)*transpose(W)*Y
%X=[temp;[biasX biasY biasZ]]
disp('where X=[ACC11 ACC21 ACC31]')
disp('        [ACC12 ACC22 ACC32]')
disp('        [ACC13 ACC23 ACC33]')
disp('        [ACC10 ACC20 ACC30]')
disp(' ');
disp(' ');
disp('calibrated=[raw,1]*A')
disp(' ');



clear W Y Y1 Y2 Y3 Y4 Y5 Y6 gyroS1 magS1 temp tempADCS1 tempS1 yawEncoder yawVelocity i
end

