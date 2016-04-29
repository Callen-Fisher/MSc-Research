%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [costFunction,costFunctionPosition,storedStates,storedPositions]=TRIAD(gyro1,gyro2,acc2,acc3,mag2,mag3,magVecX,magVecY,magVecZ,g,Q,I,covP,val1,val2,val3,val4,sampleTime,camAngles,camData,l1,l2,states,val5,val6,val7,val8)
startTime=20;
                                                                                    


R=diag([val1,val2,val3,val4]);
Q=diag([val5,val6,val7,val8]);

for i=startTime:1:length(acc2.signals.values(:,1)) 
    
    %TRIAD

g_i = [0;0;1];  % Gravitational vector in INERTIAL frame

m_i = [magVecX; magVecY; magVecZ];  % Magnetic vector in INERTIAL frame

g_b = [(acc2.signals.values(i,1)) ((acc2.signals.values(i,2))-1.08)*0 (acc2.signals.values(i,3))]'/9.81;%(acc2.signals.values(i,2))

m_b = [mag2.signals.values(i,1) mag2.signals.values(i,2) mag2.signals.values(i,3)]';

%   TRIAD algorithm

i_r = g_i/norm(g_i);
j_r = (cross(i_r, m_i))/norm(cross(i_r, m_i));
k_r = cross(i_r, j_r);

i_b = g_b/norm(g_b);
j_b = (cross(i_b, m_b))/norm(cross(i_b, m_b));
k_b = cross(i_b, j_b);

C_bn = [i_b j_b k_b]*([i_r j_r k_r]');


th1 = -1*asin(C_bn(1,3));
%ph1 = atan2(C_bn(1,2),C_bn(1,1));  % This one DOES NOT break at 90 degrees yaw
%ph1 = acos(C_bn(1,1)/cos(th1));
ph1=asin(C_bn(1,2)/cos(th1));

% if(ph1>pi/2)
%     ph1=0;
% end
% if(ph1<-pi/2)
%     ph1=0;
% end
% 
% while(ph1>pi)
%     ph1=ph1-2*pi;
% end
% while(ph1<-pi)
%     ph1=ph1+2*pi;
% end


g_b = [(acc3.signals.values(i,1)) ((acc3.signals.values(i,2))-0.6)*0 (acc3.signals.values(i,3))]'/9.81;%(acc3.signals.values(i,2))

m_b = [mag3.signals.values(i,1) mag3.signals.values(i,2) mag3.signals.values(i,3)]';

%   TRIAD algorithm

i_r = g_i/norm(g_i);
j_r = (cross(i_r, m_i))/norm(cross(i_r, m_i));
k_r = cross(i_r, j_r);

i_b = g_b/norm(g_b);
j_b = (cross(i_b, m_b))/norm(cross(i_b, m_b));
k_b = cross(i_b, j_b);

C_bn = [i_b j_b k_b]*([i_r j_r k_r]');


th2 = -1*asin(C_bn(1,3));
%ph2 = atan2(C_bn(1,2),C_bn(1,1));  % This one DOES NOT break at 90 degrees yaw
ph2=asin(C_bn(1,2)/cos(th2));

% if(ph2>pi/2)
%     ph2=0;
% end
% if(ph2<-pi/2)
%     ph2=0;
% end



% R_0_1=[[cos(states(2))*cos(states(1)) sin(states(2))*cos(states(1))   -sin(states(1))];...
%         [-sin(states(2))          cos(states(2))            0];...
%         [cos(states(2))*sin(states(1))  sin(states(2))*sin(states(1))   cos(states(1))]]';%the rotation matrix
% R_0_2=[[cos(states(4))*cos(states(3)) sin(states(4))*cos(states(3))   -sin(states(3))];...
%         [-sin(states(4))          cos(states(4))            0];...
%         [cos(states(4))*sin(states(3))  sin(states(4))*sin(states(3))   cos(states(3))]]';%the rotation matrix
%     temp=R_0_1*[0;0;ph1];
% temp2=R_0_2*[0;0;ph2];
    storedT(i,:)=[th1 ph1 th2 ph2]';
    Triad=[th1 ph1 th2 ph2];


    %KALMAN   %%%%%%%%%%IS IT MULTIPLIED BY COS TH??*cos(states(1))*cos(states(3))*cos(states(1))*cos(states(3))
    predictEq=[gyro1.signals.values(i,2) gyro1.signals.values(i,3) gyro2.signals.values(i,2) gyro2.signals.values(i,3)];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    Fmatrix=I+sampleTime*[[0 0 0 0];[0 -gyro1.signals.values(i,3)*sin(states(1)) 0 0];[0 0 0 0];[0 0 0 -gyro2.signals.values(i,3)*sin(states(3))]];%+sampleTime*[[gyro1.signals.values(i,2) 0 0 0];[0 gyro1.signals.values(i,3) 0 0];[0 0 gyro2.signals.values(i,2) 0];[0 0 0 gyro2.signals.values(i,3)]];%%%%%%%%%%%%%%%%%
    
    states = states+sampleTime*predictEq';
    
    covP = (Fmatrix)*(covP*(Fmatrix')) + Q;
   
    updateEq=states;
    
    Hmatrix=I;
    
    K=(covP*Hmatrix')/(Hmatrix*covP*Hmatrix'+R);
   
    z=Triad'; 
    states=states+K*(z-updateEq);
    covP=(I-K*Hmatrix)*covP;
    
    %store the data to be plotted
    storedStates(i,:)=states(:)';
    
    
    th1=states(1);
    ph1=states(2);
    th2=states(3);
    ph2=states(4);

R_0_1=[[cos(ph1)*cos(th1) sin(ph1)*cos(th1)   -sin(th1)];...
        [-sin(ph1)          cos(ph1)            0];...
        [cos(ph1)*sin(th1)  sin(ph1)*sin(th1)   cos(th1)]]';%the rotation matrix
R_0_2=[[cos(ph2)*cos(th2)  sin(ph2)*cos(th2)  -sin(th2)];...
        [-sin(ph2)           cos(ph2)           0];...
        [cos(ph2)*sin(th2)   sin(ph2)*sin(th2)  cos(th2)]]';

P1tip=R_0_1*[l1; 0; 0];
P2tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0];


storedPositions(i,:)=[P1tip',P2tip'];
end

%%%%%%%%%%%%%%%cam data stuff 
costFunction=sqrt((sum((storedStates(startTime:end,1)-camAngles(startTime:end,1)).^2)+...
                   sum((storedStates(startTime:end,2)-camAngles(startTime:end,2)).^2)+...
                   sum((storedStates(startTime:end,3)-camAngles(startTime:end,3)).^2)+...
                   sum((storedStates(startTime:end,4)-camAngles(startTime:end,4)).^2))/(length(acc2.signals.values(:,1))-startTime))*180/pi;
costFunctionPosition=sqrt((sum((storedPositions(startTime:end,1)-camData.signals.values(startTime:end,1)).^2)+...
                           sum((storedPositions(startTime:end,2)-camData.signals.values(startTime:end,2)).^2)+...
                           sum((storedPositions(startTime:end,3)-camData.signals.values(startTime:end,3)).^2)+...
                           sum((storedPositions(startTime:end,4)-camData.signals.values(startTime:end,4)).^2)+...
                           sum((storedPositions(startTime:end,5)-camData.signals.values(startTime:end,5)).^2)+...
                           sum((storedPositions(startTime:end,6)-camData.signals.values(startTime:end,6)).^2))/(length(acc2.signals.values(:,1))-startTime));          
x=[(storedStates(startTime:end,1)-camAngles(startTime:end,1));(storedStates(startTime:end,2)-camAngles(startTime:end,2));(storedStates(startTime:end,3)-camAngles(startTime:end,3));(storedStates(startTime:end,4)-camAngles(startTime:end,4))];
    
figure(5)
hist(x*180/pi,10)
clc  
% figure(3)
% subplot(2,2,1)
% plot(storedT(:,1))
% subplot(2,2,2)
% plot(storedT(:,2))
% subplot(2,2,3)
% plot(storedT(:,3))
% subplot(2,2,4)
% plot(storedT(:,4))

% var(storedT(:,1))
% var(storedT(:,2))
% var(storedT(:,3))
% var(storedT(:,4))
end