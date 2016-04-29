run('find_pose');

clear A1 A2 A3 A4 A_intrinsic_1 A_intrinsic_2 A_intrinsic_3 A_intrinsic_4 Rcw1 Rcw2 Rcw3 Rcw4 Tcw1 Tcw2 Tcw3 Tcw4 blue_blob cam1_blue_meas cam1_green_meas cam1_red_meas cam2_blue_meas cam2_green_meas cam2_red_meas
clear cam3_blue_meas cam3_green_meas cam3_red_meas cam4_blue_meas cam4_green_meas cam4_red_meas cam_param_1 cam_param_2 cam_param_3 cam_param_4 cx1 cx2 cx3 cx4 cy1 cy2 cy3 cy4 fd_client1 fd_client2 fd_client3 fd_client4
clear fx1 fx2 fx3 fx4 green_blob i measured object_points red_blob resnorm s1 s2 s3 s4 tl1 tl2 tl3 tl4 x x0 z fy1 fy2 fy3 fy4



%split the data into the relevant sensors:
for i=1:1:sampleLength
    s2Pos(:,i)=out(1:3,i);%red (X,Y,Z)
    s4Pos(:,i)=out(4:6,i);%green (X,Y,Z)
    s3Pos(:,i)=out(7:9,i);%blue (X,Y,Z)
end


%subtract values so that its at the "origin"
Rx=mean(s2Pos(1,1:20)); 
Ry=mean(s2Pos(2,1:20));
Rz=mean(s2Pos(3,1:20));

Bx=mean(s3Pos(1,1:20)); 
By=mean(s3Pos(2,1:20));
Bz=mean(s3Pos(3,1:20));

Gx=mean(s4Pos(1,1:20)); 
Gy=mean(s4Pos(2,1:20));
Gz=mean(s4Pos(3,1:20));

theta=atan((By-Ry)/(Bx-Rx));
theta=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=Rx;
y=Ry;
z=Rz;

originGlobal=[x-0.8,y,z]';%%%76
for i=1:1:sampleLength
    s1Pos(:,i)=[0;0;0];
    s2Pos(:,i)=s2Pos(:,i)-originGlobal;%red
    s4Pos(:,i)=s4Pos(:,i)-originGlobal;%green 
    s3Pos(:,i)=s3Pos(:,i)-originGlobal;%blue
end

R=[[cos(theta) -sin(theta) 0];[sin(theta) cos(theta) 0];[0 0 1]]';
%Rotate it so that init condition for tip is [L 0 0]
for i=1:1:sampleLength
    s2Pos(:,i)=R*s2Pos(:,i);%red
    s4Pos(:,i)=R*s4Pos(:,i);%green 
    s3Pos(:,i)=R*s3Pos(:,i);%blue
end



t=sampleLength*sampleTimeCam;
time=0:sampleTimeCam:t-sampleTimeCam;
s1Cam.time=time';
s1Cam.dimensions=3;
s1Cam.signals.values=s1Pos';
s2Cam.time=time';
s2Cam.dimensions=3;
s2Cam.signals.values=s2Pos';
s3Cam.time=time';
s3Cam.dimensions=3;
s3Cam.signals.values=s3Pos';
s4Cam.time=time';
s4Cam.dimensions=3;
s4Cam.signals.values=s4Pos';


% figure(1)
% plot3(s3(1,:),s3(2,:),s3(3,:),'r')
% hold on
% grid on
% plot3(s2(1,:),s2(2,:),s2(3,:),'b')
% plot3(0,0,0,'k+','lineWidth',5)
