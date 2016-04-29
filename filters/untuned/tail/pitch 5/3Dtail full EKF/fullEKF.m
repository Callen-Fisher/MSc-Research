%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [costFunction,costFunctionPosition,storedStates,storedPositions]=fullEKF(val1,val2,val3,val4,val5,val6,states,I,R,a1,a2,ma1,ma2,sampleTime,g1,g2,covP,camAngles,l1,l2,camData)
Q=diag([val1,val1,val2,val2,val3,val4,val5,val6,0,0,0,0]);
startTime=20;
for i=startTime:1:length(a1.signals.values(:,1)) 
    predictEq=predictFunctionFullEKF(states(1),states(3),states(2),states(4),states(6),states(8),states(5),states(7),states(10),states(12),sign(states(6)),sign(states(8)),sign(states(5)),sign(states(7)),states(9),states(11));
   
    Fmatrix=FfunctionFullEKF(states(1),states(3),states(2),states(4),states(6),states(8),states(5),states(7),states(10),states(12),sign(states(6)),sign(states(8)),sign(states(5)),sign(states(7)),states(9),states(11));
    Fmatrix=I+sampleTime*Fmatrix;
    
    states = states+sampleTime*predictEq';
    
    covP = (Fmatrix)*(covP*(Fmatrix')) + Q;
   
    updateEq=updateFunctionFullEKF(predictEq(6),predictEq(8),predictEq(5),predictEq(7),states(6),states(8),states(5),states(7),states(10),states(12),states(9),states(11));
    
    Hmatrix=HfunctionFullEKF(predictEq(6),predictEq(8),predictEq(5),predictEq(7),states(6),states(8),states(5),states(7),states(10),states(12),states(9),states(11));
    
    K=(covP*Hmatrix')/(Hmatrix*covP*Hmatrix'+R);
   
    z=[a1.signals.values(i,1) a1.signals.values(i,3) g1.signals.values(i,2) g1.signals.values(i,3) ma1.signals.values(i,1) ma1.signals.values(i,2) ma1.signals.values(i,3) a2.signals.values(i,1) a2.signals.values(i,3) g2.signals.values(i,2) g2.signals.values(i,3) ma2.signals.values(i,1) ma2.signals.values(i,2) ma2.signals.values(i,3)]'; 
    states=states+K*(z-updateEq');
    covP=(I-K*Hmatrix)*covP;
    
    %store the data to be plotted
    storedStates(i,:)=states(:)';
    
    
    th1=states(9);
ph1=states(10);
th2=states(11);
ph2=states(12);

R_0_1=[[cos(ph1) -sin(ph1) 0];[sin(ph1) cos(ph1) 0];[0 0 1]]*[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
R_0_2=[[cos(ph2) -sin(ph2) 0];[sin(ph2) cos(ph2) 0];[0 0 1]]*[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];

P1tip=R_0_1*[l1; 0; 0];
P2tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0];


storedPositions(i,:)=[P1tip',P2tip'];
end

%%%%%%%%%%%%%%%cam data stuff 
costFunction=sqrt((sum((storedStates(startTime:end,9)-camAngles(startTime:end,1)).^2)+...
                   sum((storedStates(startTime:end,10)-camAngles(startTime:end,2)).^2)+...
                   sum((storedStates(startTime:end,11)-camAngles(startTime:end,3)).^2)+...
                   sum((storedStates(startTime:end,12)-camAngles(startTime:end,4)).^2))/(length(a1.signals.values(:,1))-startTime))*180/pi;
      costFunctionPosition=sqrt((sum((storedPositions(startTime:end,1)-camData.signals.values(startTime:end,1)).^2)+...
                           sum((storedPositions(startTime:end,2)-camData.signals.values(startTime:end,2)).^2)+...
                           sum((storedPositions(startTime:end,3)-camData.signals.values(startTime:end,3)).^2)+...
                           sum((storedPositions(startTime:end,4)-camData.signals.values(startTime:end,4)).^2)+...
                           sum((storedPositions(startTime:end,5)-camData.signals.values(startTime:end,5)).^2)+...
                           sum((storedPositions(startTime:end,6)-camData.signals.values(startTime:end,6)).^2))/(length(a1.signals.values(:,1))-startTime));          
end