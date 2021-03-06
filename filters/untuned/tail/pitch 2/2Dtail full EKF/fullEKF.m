%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [costFunction,costFunctionPosition,storedStates,storedPositions]=fullEKF(val1,val2,val3,val4,states,I,R,a1,a2,sampleTime,g1,g2,covP,camAngles,ma1,ma2,l1,l2,camData)
Q=diag([val1,val2,val3,val4,0,0]);
startTime=20;
for i=startTime:1:length(a1.signals.values(:,1)) 
    predictEq=predictFunctionFullEKF(states(1),states(2),states(3),states(4),sign(states(3)),sign(states(4)),states(5),states(6));
   
    Fmatrix=FfunctionFullEKF(states(1),states(2),states(3),states(4),sign(states(3)),sign(states(4)),states(5),states(6));
    Fmatrix=I+sampleTime*Fmatrix;
    
    states = states+sampleTime*predictEq';
    
    covP = (Fmatrix)*(covP*(Fmatrix')) + Q;
   
    updateEq=updateFunctionFullEKF(predictEq(3),predictEq(4),states(3),states(4),states(5),states(6));
    
    Hmatrix=HfunctionFullEKF(predictEq(3),predictEq(4),states(3),states(4),states(5),states(6));
    
    K=(covP*Hmatrix')/(Hmatrix*covP*Hmatrix'+R);
   
    z=[a1.signals.values(i,1) a1.signals.values(i,3) g1.signals.values(i,2) ma1.signals.values(i,1) ma1.signals.values(i,3) a2.signals.values(i,1) a2.signals.values(i,3) g2.signals.values(i,2) ma2.signals.values(i,1) ma2.signals.values(i,3) ]'; 
    states=states+K*(z-updateEq');
    covP=(I-K*Hmatrix)*covP;
    
    %store the data to be plotted
    storedStates(i,:)=states(:)';
    
    
    th1=states(5);
th2=states(6);
R_0_1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
R_0_2=[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];

P1tip=R_0_1*[l1; 0; 0];
P2tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0];
    
    storedPositions(i,:)=[P1tip',P2tip'];
end

%%%%%%%%%%%%%%%cam data stuff 
costFunction=sqrt((sum((storedStates(startTime:end,5)-camAngles(startTime:end,1)).^2)+...
                   sum((storedStates(startTime:end,6)-camAngles(startTime:end,3)).^2))/(length(a1.signals.values(:,1))-startTime))*180/pi;
    costFunctionPosition=sqrt((sum((storedPositions(startTime:end,1)-camData.signals.values(startTime:end,1)).^2)+...
                           sum((storedPositions(startTime:end,3)-camData.signals.values(startTime:end,3)).^2)+...
                           sum((storedPositions(startTime:end,4)-camData.signals.values(startTime:end,4)).^2)+...
                           sum((storedPositions(startTime:end,6)-camData.signals.values(startTime:end,6)).^2))/(length(a1.signals.values(:,1))-startTime));   
            
end