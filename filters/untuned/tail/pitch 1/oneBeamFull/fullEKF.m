%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [costFunction,costFunctionPosition,storedStates,storedPositions]=fullEKF(val1,val3,val4,states,I,R,a1,a2,ma1,ma2,sampleTime,g1,g2,covP,camAngles,l1,l2,camData)
Q=diag([val1,val1,val3,val4,0,0]);
startTime=20;
for i=startTime:1:length(a1.signals.values(:,1)) 
    predictEq=predictFunctionFullEKF(states(1),states(2),states(4),states(3),states(6),sign(states(4)),sign(states(3)),states(5));%%
   
    Fmatrix=FfunctionFullEKF(states(1),states(2),states(4),states(3),states(6),sign(states(4)),sign(states(3)),states(5));%%
    Fmatrix=I+sampleTime*Fmatrix;
    
    states = states+sampleTime*predictEq';
    
    covP = (Fmatrix)*(covP*(Fmatrix')) + Q;
   
    updateEq=updateFunctionFullEKF(predictEq(2),predictEq(1),states(4),states(3),states(6),states(5));%%
    
    Hmatrix=HfunctionFullEKF(predictEq(2),predictEq(1),states(4),states(3),states(6),states(5));%%
    
    K=(covP*Hmatrix')/(Hmatrix*covP*Hmatrix'+R);
   
    z=[a2.signals.values(i,1) a2.signals.values(i,3) g2.signals.values(i,2) g2.signals.values(i,3) ma2.signals.values(i,1) ma2.signals.values(i,2) ma2.signals.values(i,3)]'; 
    states=states+K*(z-updateEq');
    covP=(I-K*Hmatrix)*covP;
    
    %store the data to be plotted
    storedStates(i,:)=states(:)';
    
    
    th=states(5);
ph=states(6);


R_0_1=[[cos(ph) -sin(ph) 0];[sin(ph) cos(ph) 0];[0 0 1]]*[[cos(th) 0 sin(th)];[0 1 0];[-sin(th) 0 cos(th)]];


P1tip=R_0_1*[l1+l2; 0; 0];


storedPositions(i,:)=[P1tip'];
end

%%%%%%%%%%%%%%%cam data stuff 
costFunction=0;
      costFunctionPosition=sqrt((sum((storedPositions(startTime:end,1)-camData.signals.values(startTime:end,4)).^2)+...
                           sum((storedPositions(startTime:end,2)-camData.signals.values(startTime:end,5)).^2)+...
                           sum((storedPositions(startTime:end,3)-camData.signals.values(startTime:end,6)).^2))/(length(a1.signals.values(:,1))-startTime));          
end