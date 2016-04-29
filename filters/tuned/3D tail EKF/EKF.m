%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [costFunction,costFunctionPosition,storedStates,storedPositions]=EKF(val3,val4,val5,val6,states,I,R,a1,a2,ma1,ma2,sampleTime,g1,g2,covP,camAngles,l1,l2,camData)
Q=diag([val3,val4,val5,val6,0,0,0,0]);
startTime=20;
for i=startTime:1:length(a1.signals.values(:,1)) 
    predictEq=predictFunctionEKF(states(2),states(4),states(1),states(3),states(6),states(8),states(5),states(7));
   
    Fmatrix=FfunctionEKF(states(2),states(4),states(1),states(3),states(6),states(8),states(5),states(7));
    Fmatrix=I+sampleTime*Fmatrix;
    
    states = states+sampleTime*predictEq';
    
    covP = (Fmatrix)*(covP*(Fmatrix')) + Q;
   
    updateEq=updateFunctionEKF(predictEq(2),predictEq(4),predictEq(1),predictEq(3),states(2),states(4),states(1),states(3),states(6),states(8),states(5),states(7));
    
    Hmatrix=HfunctionEKF(predictEq(2),predictEq(4),predictEq(1),predictEq(3),states(2),states(4),states(1),states(3),states(6),states(8),states(5),states(7));
    
    K=(covP*Hmatrix')/(Hmatrix*covP*Hmatrix'+R);
   
    z=[a1.signals.values(i,1) a1.signals.values(i,3) g1.signals.values(i,2) g1.signals.values(i,3) ma1.signals.values(i,1) ma1.signals.values(i,2) ma1.signals.values(i,3) a2.signals.values(i,1) a2.signals.values(i,3) g2.signals.values(i,2) g2.signals.values(i,3) ma2.signals.values(i,1) ma2.signals.values(i,2) ma2.signals.values(i,3)]'; 
    states=states+K*(z-updateEq');
    covP=(I-K*Hmatrix)*covP;
    
    storedStates(i,:)=states(:)';
    
    storedUpdate(i,:)=updateEq(:);
    
    
    
th1=states(5);
ph1=states(6);
th2=states(7);
ph2=states(8);

R_0_1=[[cos(ph1) -sin(ph1) 0];[sin(ph1) cos(ph1) 0];[0 0 1]]*[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
R_0_2=[[cos(ph2) -sin(ph2) 0];[sin(ph2) cos(ph2) 0];[0 0 1]]*[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];

P1tip=R_0_1*[l1; 0; 0];
P2tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0];


storedPositions(i,:)=[P1tip',P2tip'];
end

%%%%%%%%%%%%%%%cam data stuff 
costFunction=sqrt((sum((storedStates(startTime:end,5)-camAngles(startTime:end,1)).^2)+...
                   sum((storedStates(startTime:end,6)-camAngles(startTime:end,2)).^2)+...
                   sum((storedStates(startTime:end,7)-camAngles(startTime:end,3)).^2)+...
                   sum((storedStates(startTime:end,8)-camAngles(startTime:end,4)).^2))/(length(a1.signals.values(:,1))-startTime))*180/pi;
costFunctionPosition=sqrt((sum((storedPositions(startTime:end,1)-camData.signals.values(startTime:end,1)).^2)+...
                           sum((storedPositions(startTime:end,2)-camData.signals.values(startTime:end,2)).^2)+...
                           sum((storedPositions(startTime:end,3)-camData.signals.values(startTime:end,3)).^2)+...
                           sum((storedPositions(startTime:end,4)-camData.signals.values(startTime:end,4)).^2)+...
                           sum((storedPositions(startTime:end,5)-camData.signals.values(startTime:end,5)).^2)+...
                           sum((storedPositions(startTime:end,6)-camData.signals.values(startTime:end,6)).^2))/(length(a1.signals.values(:,1))-startTime));                

                       
                       
x=[(storedStates(startTime:end,5)-camAngles(startTime:end,1));(storedStates(startTime:end,6)-camAngles(startTime:end,2));(storedStates(startTime:end,7)-camAngles(startTime:end,3));(storedStates(startTime:end,8)-camAngles(startTime:end,4))];
    
figure(5)
hist(x*180/pi,10)
clc
% figure(5)
% plot(storedUpdate(:,5));
% hold on
% plot(ma1.signals.values(:,1));
% figure(6)
% plot(storedUpdate(:,6));
% hold on
% plot(ma1.signals.values(:,2));
% figure(7)
% plot(storedUpdate(:,7));
% hold on
% plot(ma1.signals.values(:,3));
end
