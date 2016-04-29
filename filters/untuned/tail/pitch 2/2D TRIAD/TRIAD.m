%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [costFunction,costFunctionPosition,storedStates,storedPositions]=TRIAD(gyro1,gyro2,acc2,acc3,mag2,mag3,magVecX,magVecY,magVecZ,g,Q,I,covP,val1,val2,sampleTime,camAngles,camData,l1,l2,states,val3,val4)
startTime=20;

R=diag([val1,val2]);
Q=diag([val3,val4]);
for i=startTime:1:length(acc2.signals.values(:,1)) 
    
    %TRIAD

g_i = [0;0;g];  % Gravitational vector in INERTIAL frame

m_i = [magVecX magVecY magVecZ]';  % Magnetic vector in INERTIAL frame

g_b = [(acc2.signals.values(i,1)) (acc2.signals.values(i,2)) (acc2.signals.values(i,3))]';

m_b = [mag2.signals.values(i,1) mag2.signals.values(i,2) mag2.signals.values(i,3)]';

%   TRIAD algorithm

i_r = g_i/norm(g_i);
j_r = (cross(i_r, m_i))/norm(cross(i_r, m_i));
k_r = cross(i_r, j_r);

i_b = g_b/norm(g_b);
j_b = (cross(i_b, m_b))/norm(cross(i_b, m_b));
k_b = cross(i_b, j_b);

C_bn = [i_b j_b k_b]*([i_r j_r k_r]');

phi = atan(C_bn(2,3)/C_bn(3,3));
theta = -1*asin(C_bn(1,3));
psi = atan2(C_bn(1,2),C_bn(1,1));  % This one DOES NOT break at 90 degrees yaw

th1=theta;

g_b = [(acc3.signals.values(i,1)) (acc3.signals.values(i,2)) (acc3.signals.values(i,3))]';

m_b = [mag3.signals.values(i,1) mag3.signals.values(i,2) mag3.signals.values(i,3)]';

%   TRIAD algorithm

i_r = g_i/norm(g_i);
j_r = (cross(i_r, m_i))/norm(cross(i_r, m_i));
k_r = cross(i_r, j_r);

i_b = g_b/norm(g_b);
j_b = (cross(i_b, m_b))/norm(cross(i_b, m_b));
k_b = cross(i_b, j_b);

C_bn = [i_b j_b k_b]*([i_r j_r k_r]');

phi = atan(C_bn(2,3)/C_bn(3,3));
theta = -1*asin(C_bn(1,3));
psi = atan2(C_bn(1,2),C_bn(1,1));  % This one DOES NOT break at 90 degrees yaw

th2=theta;


% R1wrt0=[[cos(0)*cos(th1) sin(0)*cos(th1)   -sin(th1)];...
%         [-sin(0)          cos(0)            0];...
%         [cos(0)*sin(th1)  sin(0)*sin(th1)   cos(th1)]]';%the rotation matrix
% R2wrt1=[[cos(0)*cos(th2)  sin(0)*cos(th2)  -sin(th2)];...
%         [-sin(0)           cos(0)           0];...
%         [cos(0)*sin(th2)   sin(0)*sin(th2)  cos(th2)]]';
% R2wrt0=R1wrt0*R2wrt1;
% 
% temp=R1wrt0*[0;th1;0];
% th1=temp(2);
% temp=R2wrt0*[0;th2;0];
% th2=temp(2);

    
    Triad=[th1 th2];


    %KALMAN
    predictEq=[gyro1.signals.values(i,2) gyro2.signals.values(i,2)];
   
    Fmatrix=I;%+sampleTime*[[gyro1.signals.values(i,2) 0];[0 gyro2.signals.values(i,2)]];
    
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
    th2=states(2);

R_0_1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
R_0_2=[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];

P1tip=R_0_1*[l1; 0; 0];
P2tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0];


storedPositions(i,:)=[P1tip',P2tip'];
end

%%%%%%%%%%%%%%%cam data stuff 
costFunction=sqrt((sum((storedStates(startTime:end,1)*180/pi-camAngles(startTime:end,1)*180/pi).^2)+...
                   sum((storedStates(startTime:end,2)*180/pi-camAngles(startTime:end,3)*180/pi).^2))/(length(acc2.signals.values(:,1))-startTime));
costFunctionPosition=sqrt((sum((storedPositions(startTime:end,1)-camData.signals.values(startTime:end,1)).^2)+...
                           sum((storedPositions(startTime:end,3)-camData.signals.values(startTime:end,3)).^2)+...
                           sum((storedPositions(startTime:end,4)-camData.signals.values(startTime:end,4)).^2)+...
                           sum((storedPositions(startTime:end,6)-camData.signals.values(startTime:end,6)).^2))/(length(acc2.signals.values(:,1))-startTime));          
end