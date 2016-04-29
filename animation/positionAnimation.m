clc
clear

%TODO
%legend

%this must be removed and received in the function
% fileFilter='positionSpine';
% fileCamera='positionCamera';
fileFilter='pos';%one sensor, one beam
fileCamera='cam';
vid2Dname='twoD.avi';
vid3Dname='threeD.avi';
%fileFilter='positionTail';
%fileCamera='camTail';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fileFilter);
load(fileCamera);

%remove the data from its structure (easier to work with)
Pfilter=storedPositions; %data length X positions
Pcamera=camData.signals.values(:,:);

clear camData storedPositions fileCamera fileFilter

if(length(Pcamera(1,:))==12)
    %spine (12 positions)
    type=1;
elseif(length(Pfilter(1,:))==6)
    %tail (6 positions (X,Y,Z of the 2 blobs (base it at the origin))
    type=0;
else
    %one beam one sensor tail
    type=2;
end

%3D animation
l1=0.4;
l2=0.3;

writerObj = VideoWriter(vid3Dname);%name of the video file 
open(writerObj);%opens the object
warning('off');
figure(1)
for i=1:1:length(Pcamera(:,1))
    clf;
    axis([-0.5 2 -0.5 0.5 -2 2]);
    hold on
    xlabel('X axis');
    ylabel('Y axis');
    zlabel('Z axis');
    if(type==1)
        title('3D image of tail and spine flick')
        l3=sqrt(Pcamera(i,4)^2+Pcamera(i,5)^2+Pcamera(i,6)^2);
        th1=asin(Pcamera(i,6)/l3);
        ph=asin(Pcamera(i,5)/l3);
        th2=acos((l1^2+l3^3-l2^2)/(2*l1*l3));
        th=-(th1+th2);
        Rpitch=[[cos(th) 0 sin(th)];[0 1 0];[-sin(th) 0 cos(th)]];
        Ryaw=[[cos(ph) -sin(ph) 0];[sin(ph) cos(ph) 0];[0 0 1]];
        R=Ryaw*Rpitch;
        Px=R*[l1;0;0];
        
        l3=sqrt(Pfilter(i,4)^2+Pfilter(i,5)^2+Pfilter(i,6)^2);
        th1=asin(Pfilter(i,6)/l3);
        ph=asin(Pfilter(i,5)/l3);
        th2=acos((l1^2+l3^3-l2^2)/(2*l1*l3));
        th=-(th1+th2);
        Rpitch=[[cos(th) 0 sin(th)];[0 1 0];[-sin(th) 0 cos(th)]];
        Ryaw=[[cos(ph) -sin(ph) 0];[sin(ph) cos(ph) 0];[0 0 1]];
        R=Ryaw*Rpitch;
        Px2=R*[l1;0;0];
        
        %plot camera
        plot3([Pcamera(i,1),Px(1),Pcamera(i,4),Pcamera(i,7),Pcamera(i,10)],[Pcamera(i,2),Px(2),Pcamera(i,5),Pcamera(i,8),Pcamera(i,11)],[Pcamera(i,3),Px(3),Pcamera(i,6),Pcamera(i,9),Pcamera(i,12)]);
        %plot filter
        plot3([Pfilter(i,1),Px2(1),Pfilter(i,4),Pfilter(i,7),Pfilter(i,10)],[Pfilter(i,2),Px2(2),Pfilter(i,5),Pfilter(i,8),Pfilter(i,11)],[Pfilter(i,3),Px2(3),Pfilter(i,6),Pfilter(i,9),Pfilter(i,12)],'r');
        legend('Camera system','Filter system');
    elseif(type==0)
        title('3D image of tail flick')
        %plot camera 
        plot3([0,Pcamera(i,1),Pcamera(i,4)],[0,Pcamera(i,2),Pcamera(i,5)],[0,Pcamera(i,3),Pcamera(i,6)]);
        %plot filter
        plot3([0,Pfilter(i,1),Pfilter(i,4)],[0,Pfilter(i,2),Pfilter(i,5)],[0,Pfilter(i,3),Pfilter(i,6)],'r');
        legend('Camera system','Filter system');
    elseif(type==2)
        plot3([0,Pcamera(i,1),Pcamera(i,4)],[0,Pcamera(i,2),Pcamera(i,5)],[0,Pcamera(i,3),Pcamera(i,6)]);
        plot3([0,Pfilter(i,1)],[0,Pfilter(i,2)],[0,Pfilter(i,3)],'r');
    end
    plot3(2,2,2)
    plot3(-0.5,-2,-2)
    frame = getframe(gcf);%make the movie
    writeVideo(writerObj,frame);
    pause(0.01);
end
close(writerObj);

%2D animation
writerObj = VideoWriter(vid2Dname);%name of the video file 
open(writerObj);%opens the object
figure(2)
for i=1:1:length(Pcamera(:,1))
    clf;
    axis([-0.5 2 -2 2])
    hold on
    xlabel('X axis');
    ylabel('Y axis');
    if(type==1)
        title('2D image of tail and spine flick')
        %plot camera
        l3=sqrt(Pcamera(i,4)^2+Pcamera(i,5)^2+Pcamera(i,6)^2);
        th1=asin(Pcamera(i,6)/l3);
        ph=asin(Pcamera(i,5)/l3);
        th2=acos((l1^2+l3^3-l2^2)/(2*l1*l3));
        th=-(th1+th2);
        Rpitch=[[cos(th) 0 sin(th)];[0 1 0];[-sin(th) 0 cos(th)]];
        Ryaw=[[cos(ph) -sin(ph) 0];[sin(ph) cos(ph) 0];[0 0 1]];
        R=Ryaw*Rpitch;
        Px=R*[l1;0;0];
        
        l3=sqrt(Pfilter(i,4)^2+Pfilter(i,5)^2+Pfilter(i,6)^2);
        th1=asin(Pfilter(i,6)/l3);
        ph=asin(Pfilter(i,5)/l3);
        th2=acos((l1^2+l3^3-l2^2)/(2*l1*l3));
        th=-(th1+th2);
        Rpitch=[[cos(th) 0 sin(th)];[0 1 0];[-sin(th) 0 cos(th)]];
        Ryaw=[[cos(ph) -sin(ph) 0];[sin(ph) cos(ph) 0];[0 0 1]];
        R=Ryaw*Rpitch;
        Px2=R*[l1;0;0];
        
        plot([Pcamera(i,1),Px(1),Pcamera(i,4),Pcamera(i,7),Pcamera(i,10)],[Pcamera(i,3),Px(3),Pcamera(i,6),Pcamera(i,9),Pcamera(i,12)]);
        %plot filter
        plot([Pfilter(i,1),Px2(1),Pfilter(i,4),Pfilter(i,7),Pfilter(i,10)],[Pfilter(i,3),Px(3),Pfilter(i,6),Pfilter(i,9),Pfilter(i,12)],'r');
        legend('Camera system','Filter system');
    elseif(type==0)
        title('2D image of tail flick')
        %plot camera 
        plot([0,Pcamera(i,1),Pcamera(i,4)],[0,Pcamera(i,3),Pcamera(i,6)]);
        %plot filter
        plot([0,Pfilter(i,1),Pfilter(i,4)],[0,Pfilter(i,3),Pfilter(i,6)],'r');
        legend('Camera system','Filter system');
    else
        plot([0,Pcamera(i,1),Pcamera(i,4)],[0,Pcamera(i,3),Pcamera(i,6)]);
        plot([0,Pfilter(i,1)],[0,Pfilter(i,3)],'r');
    end
    plot(2,2)
    plot(-0.5,-2)
    frame = getframe(gcf);%make the movie
    writeVideo(writerObj,frame);
    pause(0.01);
end
close(writerObj);

clear all
disp('FINISHED')