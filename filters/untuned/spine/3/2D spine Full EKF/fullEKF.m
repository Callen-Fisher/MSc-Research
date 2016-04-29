%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [costFunction,costFunctionPosition,storedStates,storedPositions]=fullEKF(val1,val2,val3,val4,val5,val6,states,I,R,a1,a2,a3,a4,sampleTime,g1,g2,g3,g4,covP,camAngles,ma1,ma2,ma3,ma4,l1,l2,l3,l4,camData)
Q=diag([val1,val2,val3,val4,val5,val6,0,0,0,0]);
startTime=20;
for i=startTime:1:length(a1.signals.values(:,1)) 
    predictEq=predictFunctionFullEKF(states(1),states(2),states(3),states(4),states(5),states(6),sign(5),sign(6),states(7),states(8),states(9),states(10));
   
    Fmatrix=FfunctionFullEKF(states(1),states(2),states(3),states(4),states(5),states(6),sign(5),sign(6),states(7),states(8),states(9),states(10));
    Fmatrix=I+sampleTime*Fmatrix;
    
    states = states+sampleTime*predictEq';
    
    covP = (Fmatrix)*(covP*(Fmatrix')) + Q;
   
    updateEq=updateFunctionFullEKF(predictEq(1),predictEq(2),predictEq(3),predictEq(4),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10));
    
    Hmatrix=HfunctionFullEKF(predictEq(1),predictEq(2),predictEq(3),predictEq(4),states(3),states(4),states(5),states(6),states(7),states(8),states(9),states(10));
    
    K=(covP*Hmatrix')/(Hmatrix*covP*Hmatrix'+R);
   
    z=[a1.signals.values(i,1) a1.signals.values(i,3) g1.signals.values(i,2) ma1.signals.values(i,1) ma1.signals.values(i,3)...
       a2.signals.values(i,1) a2.signals.values(i,3) g2.signals.values(i,2) ma2.signals.values(i,1) ma2.signals.values(i,3)...
       a3.signals.values(i,1) a3.signals.values(i,3) g3.signals.values(i,2) ma3.signals.values(i,1) ma3.signals.values(i,3)...
       a4.signals.values(i,1) a4.signals.values(i,3) g4.signals.values(i,2) ma4.signals.values(i,1) ma4.signals.values(i,3)]'; 
    states=states+K*(z-updateEq');
    covP=(I-K*Hmatrix)*covP;
    
    storedStates(i,:)=states(:)';
    
    
    
    
th1=states(7);
th2=states(8);
th3=states(9);
th4=states(10);
R_0_1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
R_0_2=[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];
R_0_3=[[cos(th3) 0 sin(th3)];[0 1 0];[-sin(th3) 0 cos(th3)]];
R_0_4=[[cos(th4) 0 sin(th4)];[0 1 0];[-sin(th4) 0 cos(th4)]];

P1base=[0;0;0];
P2tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0];
P3tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[l3;0;0];
P4tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[l3;0;0]+R_0_4*[l4;0;0];
    
    
    storedPositions(i,:)=[P1base',P2tip',P3tip',P4tip'];
end

%%%%%%%%%%%%%%%cam data stuff 
costFunction=sqrt((sum((storedStates(startTime:end,7)-camAngles(startTime:end,1)).^2)+...
                   sum((storedStates(startTime:end,8)-camAngles(startTime:end,2)).^2)+...
                   sum((storedStates(startTime:end,9)-camAngles(startTime:end,3)).^2)+...
                   sum((storedStates(startTime:end,10)-camAngles(startTime:end,4)).^2))/(length(a1.signals.values(:,1))-startTime))*180/pi;
costFunctionPosition=sqrt((sum((storedPositions(startTime:end,1)-camData.signals.values(startTime:end,1)).^2)+...
                           sum((storedPositions(startTime:end,3)-camData.signals.values(startTime:end,3)).^2)+...
                           sum((storedPositions(startTime:end,4)-camData.signals.values(startTime:end,4)).^2)+...
                           sum((storedPositions(startTime:end,6)-camData.signals.values(startTime:end,6)).^2)+...
                           sum((storedPositions(startTime:end,7)-camData.signals.values(startTime:end,7)).^2)+...
                           sum((storedPositions(startTime:end,9)-camData.signals.values(startTime:end,9)).^2)+...
                           sum((storedPositions(startTime:end,10)-camData.signals.values(startTime:end,10)).^2)+...
                           sum((storedPositions(startTime:end,12)-camData.signals.values(startTime:end,12)).^2))/(length(a1.signals.values(:,1))-startTime));   
        
end
