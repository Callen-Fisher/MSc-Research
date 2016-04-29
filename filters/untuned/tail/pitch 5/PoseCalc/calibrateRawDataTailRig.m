%fileName='testSingleHandFlick.dat';
run('importDataNsensors');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the motor torque
TorqueMotor=[];
torqueScalingFactor=100;
noCurrentValue=100;

for i=1:1:length(tempS4(:,1))%%%%THIS MUST CHANGE
    if(tempS4(i,1)<1)
        tempS4(i,1)=noCurrentValue;
    end
    TorqueMotor(i)=(tempS4(i,1)-noCurrentValue)*torqueScalingFactor;
end


%clear gyroS4 magS4 accS4 tempS4 tempADCS4 torqueScalingFactor noCurrentValue fileName
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the Acc for S1
%calibrated=[raw,1]*A

A=[ [0.9576   -0.0014   -0.0058];
    [0.0084    0.9446    0.0063];
    [-0.0049   -0.0024    0.9524];
    [0.0226    0.0105    0.1726]];

for i=1:1:length(accS1(:,1))
    calAccS1(i,:)=double([accS1(i,:),1]*A);
end

%get the Acc for S2
%calibrated=[raw,1]*A


A=[ [0.9217   -0.0237    0.0055];
    [0.0337    0.9665   -0.0120];
    [-0.0173    0.0037    0.9211];
    [0.0482   -0.0037    0.0720]];

for i=1:1:length(accS2(:,1))
    calAccS2(i,:)=double([accS2(i,:),1]*A);%*0.8369;%%%%%%%%%%%%%%%%%%%%%
end
%get the Acc for S3
%[Ax1] [ACC11 ACC12 ACC13][Ax] [ACC10]
%[Ay1]=[ACC21 ACC22 ACC23][Ay]+[ACC20]
%[Az1] [ACC31 ACC32 ACC33][Az] [ACC30]

%new=A*(raw-bias)
A=[ [0.9243    0.0029   -0.0191];
    [-0.0063    0.9514    0.0537];
    [-0.0206    0.0012    0.9107];
    [0.0430    0.0129    0.1162]];


for i=1:1:length(accS3(:,1))
    calAccS3(i,:)=double([accS3(i,:),1]*A);%*0.8056;%%%%%%%%%%%%%%%
end


%get the acc for S4
A=[ [0.1543    0.0025    0.0003];
    [-0.0006    0.1591    0.0028];
    [-0.0017    0.0010    0.1552];
    [0.0137   -0.0152    0.0957]];


for i=1:1:length(accS4(:,1))
    calAccS4(i,:)=double([accS4(i,:),1]*A);%*3.9103;%%%%%%%%%%%%%%%
end




clear accS3 accS2 accS1 accS4 A bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the Mag for S1
%calibrated measurement w = U*(v-c)');
U=[ [ 4.0534    0.4806    0.0899];
    [ 0    5.0900    0.2135];
    [0         0    3.6083]];
c=[0.1707 0.1738 -0.1053]';

calMagS1(1,:)=double((U*(magS1(1,:)-c')')');
for i=2:1:length(magS1(:,1))
    calMagS1(i,:)=double((U*(magS1(i,:)-c')')');
    if(calMagS1(i,1)^2>1.2)
        calMagS1(i,1)=calMagS1(i-1,1);
    end
    if(calMagS1(i,2)^2>1.2)
        calMagS1(i,2)=calMagS1(i-1,2);
    end
    if(calMagS1(i,3)^2>1.2)
        calMagS1(i,3)=calMagS1(i-1,3);
    end
end
%get the Mag for S2
%calibrated measurement w = U*(v-c)');
U=[ [3.535    0.2141    -0.2101];
    [0    4.5164   -0.2115];
    [0         0    3.5389]];
c=[0.0346 -0.0231 -0.2341]';

calMagS2(1,:)=double((U*(magS2(1,:)-c')')');
for i=2:1:length(magS2(:,1))
    calMagS2(i,:)=double((U*(magS2(i,:)-c')')');
    if(calMagS2(i,1)^2>1.2)
        calMagS2(i,1)=calMagS2(i-1,1);
    end
    if(calMagS2(i,2)^2>1.2)
        calMagS2(i,2)=calMagS2(i-1,2);
    end
    if(calMagS2(i,3)^2>1.2)
        calMagS2(i,3)=calMagS2(i-1,3);
    end
end
%get the Mag for S3
%calibrated measurement w = U*(v-c)');
U=[ [3.7927    0.1502   -0.0791];
    [0    4.4416    0.0277];
    [0         0    3.2403]];
c=[0.0091 -0.0239 -0.1802]';

calMagS3(1,:)=double((U*(magS3(1,:)-c')')');
for i=2:1:length(magS3(:,1))
    calMagS3(i,:)=double((U*(magS3(i,:)-c')')');
    if(calMagS3(i,1)^2>1.2)
        calMagS3(i,1)=calMagS3(i-1,1);
    end
    if(calMagS3(i,2)^2>1.2)
        calMagS3(i,2)=calMagS3(i-1,2);
    end
    if(calMagS3(i,3)^2>1.2)
        calMagS3(i,3)=calMagS3(i-1,3);
    end
end

%get the Mag for S4
%calibrated measurement w = U*(v-c)');
U=[ [3.5589    0.0788   -0.0108];
    [0    4.4166    0.0803];
    [0         0    3.1110]];
c=[0.0003 0.0058 -0.1819]';
calMagS4(1,:)=double((U*(magS4(1,:)-c')')');%;*0.2748;%%%%%%%%%%%%%%%%
for i=2:1:length(magS4(:,1))
    calMagS4(i,:)=double((U*(magS4(i,:)-c')')');%;*0.2748;%%%%%%%%%%%%%%%
%     if(calMagS4(i,1)^2>1.2)
%         calMagS4(i,1)=calMagS4(i-1,1);
%     end
%     if(calMagS4(i,2)^2>1.2)
%         calMagS4(i,2)=calMagS4(i-1,2);
%     end
%     if(calMagS4(i,3)^2>1.2)
%         calMagS4(i,3)=calMagS4(i-1,3);
%     end
end
clear magS3 magS2 magS1  magS4 U c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the gyro for S1
grad1=-0.0021;
grad2=-0.0087;
grad3=0.031;
offset1=-1.7;
offset2=+0.74;
offset3=+1.2;

for i=1:1:length(gyroS1(:,1))
    calgyroS1(i,1)=double(gyroS1(i,1)-offset1-grad1*single(tempADCS1(i)));
    calgyroS1(i,2)=double(gyroS1(i,2)-offset2-grad2*single(tempADCS1(i)));
    calgyroS1(i,3)=double(gyroS1(i,3)-offset3-grad3*single(tempADCS1(i)));
end
%get the gyro for S2
grad1=-0.041;
grad2=-0.0094;
grad3=-0.021;
offset1=-1.9;
offset2=+1.6;
offset3=+0.11;

for i=1:1:length(gyroS2(:,1))
    calgyroS2(i,1)=double(gyroS2(i,1)-offset1-grad1*single(tempADCS2(i)));
    calgyroS2(i,2)=double(gyroS2(i,2)-offset2-grad2*single(tempADCS2(i)));
    calgyroS2(i,3)=double(gyroS2(i,3)-offset3-grad3*single(tempADCS2(i)));
end
%get the gyro for S3
grad1=-0.0053;
grad2=-0.025;
grad3=-0.0045;
offset1=+0.94;
offset2=-5.3;
offset3=-1.1;

for i=1:1:length(gyroS3(:,1))
    calgyroS3(i,1)=double(gyroS3(i,1)-offset1-grad1*single(tempADCS3(i)));
    calgyroS3(i,2)=double(gyroS3(i,2)-offset2-grad2*single(tempADCS3(i)));
    calgyroS3(i,3)=double(gyroS3(i,3)-offset3-grad3*single(tempADCS3(i)));
end
%get the gyro for S4
grad1=-0.011;
grad2=-0.0017;
grad3=-0.00056;
offset1=-0.37;
offset2=+0.91;
offset3=-1;

for i=1:1:length(gyroS1(:,1))
    calgyroS4(i,1)=double(gyroS4(i,1)-offset1-grad1*single(tempADCS4(i)));
    calgyroS4(i,2)=double(gyroS4(i,2)-offset2-grad2*single(tempADCS4(i)));
    calgyroS4(i,3)=double(gyroS4(i,3)-offset3-grad3*single(tempADCS4(i)));
end


t=length(gyroS3(:,1))*0.01;

%align the axes and put it into structures:
degToRad=pi/180;
time=0:0.01:t-0.01;
DimValues=3;
%the rotation of the sensor data
gyro1.time=time';
gyro1.dimensions=3;
gyro1.signals.values=[-calgyroS1(:,1),-calgyroS1(:,2),calgyroS1(:,3)]*degToRad;
gyro2.time=time';
gyro2.dimensions=3;
gyro2.signals.values=[-calgyroS2(:,1),-calgyroS2(:,2),calgyroS2(:,3)]*degToRad;
gyro3.time=time';
gyro3.dimensions=3;
gyro3.signals.values=[-calgyroS3(:,1),-calgyroS3(:,2),calgyroS3(:,3)]*degToRad;
gyro4.time=time';
gyro4.dimensions=3;
gyro4.signals.values=[-calgyroS4(:,1),-calgyroS4(:,2),calgyroS4(:,3)]*degToRad;

acc1.time=time';
acc1.dimensions=3;
acc1.signals.values=[-calAccS1(:,2),calAccS1(:,1),calAccS1(:,3)]*9.81;
acc2.time=time';
acc2.dimensions=3;
acc2.signals.values=[-calAccS2(:,2),calAccS2(:,1),calAccS2(:,3)]*9.81;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acc3.time=time';
acc3.dimensions=3;
acc3.signals.values=[-calAccS3(:,2),calAccS3(:,1),calAccS3(:,3)]*9.81;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acc4.time=time';
acc4.dimensions=3;
acc4.signals.values=[-calAccS4(:,2),calAccS4(:,1),calAccS4(:,3)]*9.81;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mag1.time=time';
mag1.dimensions=3;
mag1.signals.values=[-calMagS1(:,2),calMagS1(:,1),calMagS1(:,3)];
mag2.time=time';
mag2.dimensions=3;
mag2.signals.values=[-calMagS2(:,2),calMagS2(:,1),calMagS2(:,3)];
mag3.time=time';
mag3.dimensions=3;
mag3.signals.values=[-calMagS3(:,2),calMagS3(:,1),calMagS3(:,3)];
mag4.time=time';
mag4.dimensions=3;
mag4.signals.values=[-calMagS4(:,2),calMagS4(:,1),calMagS4(:,3)];








% 
% angle1=10*pi/180;
% R1=[[1 0 0];[0 cos(angle1) -sin(angle1)];[0 sin(angle1) cos(angle1)]]';
% angle2=5*pi/180;
% R2=[[cos(angle2) -sin(angle2) 0];[sin(angle2) cos(angle2) 0];[0 0 1]]';
% R = R2*R1;
% for i=1:1:length(acc2.signals.values(:,1))
%     %rotate the acc vector:
%    temp=R'*acc2.signals.values(i,:)';
%    acc2.signals.values(i,1)=temp(1);
%    acc2.signals.values(i,2)=temp(2);
%    acc2.signals.values(i,3)=temp(3);
%    
%    temp=R'*acc3.signals.values(i,:)';
%    acc3.signals.values(i,1)=temp(1);
%    acc3.signals.values(i,2)=temp(2);
%    acc3.signals.values(i,3)=temp(3);
%  
%    %rotate the mag vector:
%    temp=R*mag2.signals.values(i,:)';
%    mag2.signals.values(i,1)=temp(1);
%    mag2.signals.values(i,2)=temp(2);
%    mag2.signals.values(i,3)=temp(3);
%    
%    temp=R*mag3.signals.values(i,:)';
%    mag3.signals.values(i,1)=temp(1);
%    mag3.signals.values(i,2)=temp(2);
%    mag3.signals.values(i,3)=temp(3);
%  
%    %rotate the gyro
%    temp=R*gyro2.signals.values(i,:)';
%    gyro2.signals.values(i,1)=temp(1);
%    gyro2.signals.values(i,2)=temp(2);
%    gyro2.signals.values(i,3)=temp(3);
%    
%    temp=R*gyro3.signals.values(i,:)';
%    gyro3.signals.values(i,1)=temp(1);
%    gyro3.signals.values(i,2)=temp(2);
%    gyro3.signals.values(i,3)=temp(3);
%  
% end










% angle=85*pi/180;
% R=[[cos(angle) 0 sin(angle)];[0 1 0];[-sin(angle) 0 cos(angle)]]
% %rotate back up to horizontal position
% for i=1:1:length(acc2.signals.values(:,1))
%     %rotate the acc vector:
%    temp=R*acc2.signals.values(i,:)';
%    acc2.signals.values(i,1)=temp(1);
%    acc2.signals.values(i,2)=temp(2);
%    acc2.signals.values(i,3)=temp(3);
%    
%    temp=R*acc3.signals.values(i,:)';
%    acc3.signals.values(i,1)=temp(1);
%    acc3.signals.values(i,2)=temp(2);
%    acc3.signals.values(i,3)=temp(3);
%  
%    %rotate the mag vector:
%    temp=R*mag2.signals.values(i,:)';
%    mag2.signals.values(i,1)=temp(1);
%    mag2.signals.values(i,2)=temp(2);
%    mag2.signals.values(i,3)=temp(3);
%    
%    temp=R*mag3.signals.values(i,:)';
%    mag3.signals.values(i,1)=temp(1);
%    mag3.signals.values(i,2)=temp(2);
%    mag3.signals.values(i,3)=temp(3);
%  
%    %rotate the gyro
%    temp=R*gyro2.signals.values(i,:)';
%    gyro2.signals.values(i,1)=temp(1);
%    gyro2.signals.values(i,2)=temp(2);
%    gyro2.signals.values(i,3)=temp(3);
%    
%    temp=R*gyro3.signals.values(i,:)';
%    gyro3.signals.values(i,1)=temp(1);
%    gyro3.signals.values(i,2)=temp(2);
%    gyro3.signals.values(i,3)=temp(3);
%  
% end
% 
% 
% 
% 






% %calculate the roll:
% lengthAcc2=0;
% 
% tempX=0;
% tempY=0;
% tempZ=0;
% 
% for i=1:1:10
%     lengthAcc2=sqrt(acc2.signals.values(i,1)^2+acc2.signals.values(i,2)^2+acc2.signals.values(i,3)^2);
%     tempX=acc2.signals.values(i,1)+tempX;
%     tempY=acc2.signals.values(i,2)+tempY;
%     tempZ=acc2.signals.values(i,3)+tempZ;
% end
% lengthAcc2=lengthAcc2/10;
% tempX=tempX/10;
% tempY=tempY/10;
% tempZ=tempZ/10;
% 
% rollS2=atan2(tempY,tempZ)
% rollS2=-20*pi/180
% lengthAcc3=0;
% for i=1:1:10
%     lengthAcc3=sqrt(acc3.signals.values(i,1)^2+acc3.signals.values(i,2)^2+acc3.signals.values(i,3)^2);
%     tempX=acc3.signals.values(i,1)+tempX;
%     tempY=acc3.signals.values(i,2)+tempY;
%     tempZ=acc3.signals.values(i,3)+tempZ;
% end
% lengthAcc3=lengthAcc3/10;
% tempX=tempX/10;
% tempY=tempY/10;
% tempZ=tempZ/10;
% 
% rollS3=atan2(tempY,tempZ)
% 
% R2=[[1 0 0];[0 cos(-rollS2) -sin(-rollS2)];[0 sin(-rollS2) cos(-rollS2)]];
% R3=[[1 0 0];[0 cos(-rollS3) -sin(-rollS3)];[0 sin(-rollS3) cos(-rollS3)]];
% 
% for i=1:1:length(acc2.signals.values(:,1))
%     %rotate the acc vector:
%    temp=R2*acc2.signals.values(i,:)';
%    acc2.signals.values(i,1)=temp(1);
%    acc2.signals.values(i,2)=temp(2);
%    acc2.signals.values(i,3)=temp(3);
%    
%    temp=R3*acc3.signals.values(i,:)';
%    acc3.signals.values(i,1)=temp(1);
%    acc3.signals.values(i,2)=temp(2);
%    acc3.signals.values(i,3)=temp(3);
%  
%    %rotate the mag vector:
%    temp=R2*mag2.signals.values(i,:)';
%    mag2.signals.values(i,1)=temp(1);
%    mag2.signals.values(i,2)=temp(2);
%    mag2.signals.values(i,3)=temp(3);
%    
%    temp=R3*mag3.signals.values(i,:)';
%    mag3.signals.values(i,1)=temp(1);
%    mag3.signals.values(i,2)=temp(2);
%    mag3.signals.values(i,3)=temp(3);
%  
%    %rotate the gyro
%    temp=R2*gyro2.signals.values(i,:)';
%    gyro2.signals.values(i,1)=temp(1);
%    gyro2.signals.values(i,2)=temp(2);
%    gyro2.signals.values(i,3)=temp(3);
%    
%    temp=R3*gyro3.signals.values(i,:)';
%    gyro3.signals.values(i,1)=temp(1);
%    gyro3.signals.values(i,2)=temp(2);
%    gyro3.signals.values(i,3)=temp(3);
%  
% end
% 
% R=[[cos(-angle) 0 sin(-angle)];[0 1 0];[-sin(-angle) 0 cos(-angle)]]
% %rotate back up to vertical position
% for i=1:1:length(acc2.signals.values(:,1))
%     %rotate the acc vector:
%    temp=R*acc2.signals.values(i,:)';
%    acc2.signals.values(i,1)=temp(1);
%    acc2.signals.values(i,2)=temp(2);
%    acc2.signals.values(i,3)=temp(3);
%    
%    temp=R*acc3.signals.values(i,:)';
%    acc3.signals.values(i,1)=temp(1);
%    acc3.signals.values(i,2)=temp(2);
%    acc3.signals.values(i,3)=temp(3);
%  
%    %rotate the mag vector:
%    temp=R*mag2.signals.values(i,:)';
%    mag2.signals.values(i,1)=temp(1);
%    mag2.signals.values(i,2)=temp(2);
%    mag2.signals.values(i,3)=temp(3);
%    
%    temp=R*mag3.signals.values(i,:)';
%    mag3.signals.values(i,1)=temp(1);
%    mag3.signals.values(i,2)=temp(2);
%    mag3.signals.values(i,3)=temp(3);
%  
%    %rotate the gyro
%    temp=R*gyro2.signals.values(i,:)';
%    gyro2.signals.values(i,1)=temp(1);
%    gyro2.signals.values(i,2)=temp(2);
%    gyro2.signals.values(i,3)=temp(3);
%    
%    temp=R*gyro3.signals.values(i,:)';
%    gyro3.signals.values(i,1)=temp(1);
%    gyro3.signals.values(i,2)=temp(2);
%    gyro3.signals.values(i,3)=temp(3);
%  
% end



torqueM.signals.values=TorqueMotor';
torqueM.dimensions=1;
torqueM.time=time';

clear gyroS1 gyroS2 gyroS3 gyroS4 tempADCS1 tempADCS2 tempADCS3 tempADCS4 tempS1 tempS2 tempS3 tempS4
clear grad1 grad2 grad3 offset1 offset2 offset3 i

clear accScale gyroScale magScale1 magScale2 tempADCScale tempScale