clc


for i=1:1:length(g3.signals.values(:,1))
    p1=camData.signals.values(i,1:3)';
    p2=camData.signals.values(i,4:6)';
    p3=camData.signals.values(i,7:9)';
    p4=camData.signals.values(i,10:12)';
    
    lx=sqrt((p2(1)-p1(1))^2+(p2(3)-p1(3))^2);
    temp=(-l2^2+l1^2+lx^2)/(2*l1*lx);
    if(temp^2>1.0)
        if(temp>1)
            temp=1;
        else
            temp=-1;
        end
    end
    th1_x=acos(temp);
    thx=asin((p1(3)-p2(3))/lx);
    th1=-th1_x+thx;%%%%%%%%neg
    
    rot=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
    Pmid=rot*[l1;0;0];
    th2=-atan2((p2(3)-Pmid(3)),(p2(1)-Pmid(1)))+2*pi/180;
    
    th1=th1;
    th2=th2;
    
    th3=-atan2((p3(3)-p2(3)),(p3(1)-p2(1)));
    th4=-atan2((p4(3)-p3(3)),(p4(1)-p3(1)));

    camAngles(i,:)=[th1,th2,th3,th4];
end


% figure(1)
% plot(camAngles(:,1));
% figure(2)
% plot(camAngles(:,2));