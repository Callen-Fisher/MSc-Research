clc


for i=1:1:length(g1.signals.values(:,1))
    p1=camData.signals.values(i,1:3)';
    p2=camData.signals.values(i,4:6)';
    
    th1=-atan2(p1(3),p1(1));
    th2=-atan2((p2(3)-p1(3)),(p2(1)-p1(1)));

    th1=th1+10*pi/180;
    th2=th2+10*pi/180;
    %rotate the point down!!!!
    ph1=0;
    ph2=0;
    R1wrt0T=[[cos(ph1)*cos(th1) sin(ph1)*cos(th1)   -sin(th1)];...
            [-sin(ph1)          cos(ph1)            0];...
            [cos(ph1)*sin(th1)  sin(ph1)*sin(th1)   cos(th1)]];%the rotation matrix
    R2wrt0T=[[cos(ph2)*cos(th2)  sin(ph2)*cos(th2)  -sin(th2)];...
            [-sin(ph2)           cos(ph2)           0];...
            [cos(ph2)*sin(th2)   sin(ph2)*sin(th2)  cos(th2)]];


    p1=R1wrt0T*p1;
    p2=R2wrt0T*p2;

    ph1=atan2(p1(2),p1(1));
    ph2=atan2((p2(2)-p1(2)),(p2(1)-p1(1)));
    ph2temp=R1wrt0T*[0;0;ph2];
    ph2=ph2temp(3);

    camAngles(i,:)=[th1,ph1,th2,ph2];
end