function [costFunction,storedStates,storedPositions]=EKF(val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,states,covP,R,I,sampleTime,camData,a1,a2,a3,a4,g1,g2,g3,g4,ma1,ma2,ma3,ma4,l1,l2,l3,l4)
Q_noise=diag([val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,0,0,0,0,0,0,0,0,0,0]);
startTime=20;
warning('off');
for i=startTime:1:length(a1.signals.values(:,1)) 
    Q=[0;0;0;0;0;0;0;0;0;0];
    %M=PH1,PH2,PH3,PH4,TH1,TH2,TH3,TH4
    M=MfunctionEKF(states(13),states(16),states(18),states(20),states(12),states(15),states(17),states(19));
    %C=DPH1,DPH2,DPH3,DPH4,DSI1,DSI2,DTH1,DTH2,DTH3,DTH4,PH1,PH2,PH3,PH4,TH1,TH2,TH3,TH4
    C=CfunctionEKF(states(3),states(6),states(8),states(10),states(1),states(4),states(2),states(5),states(7),states(9),...
        states(13),states(16),states(18),states(20),states(12),states(15),states(17),states(19));
    dq=[states(1); states(2); states(3); states(4); states(5); states(6); states(7); states(8); states(9); states(10)];
    G=[0;0;0;0;0;0;0;0;0;0];
    systemEquations=(M)\(-C*dq - G + Q);
    predictEq=[systemEquations(1),systemEquations(2),systemEquations(3),systemEquations(4),systemEquations(5),systemEquations(6),...
               systemEquations(7),systemEquations(8),systemEquations(9),systemEquations(10),states(1),states(2),states(3),states(4),...
               states(5),states(6),states(7),states(8),states(9),states(10)];
    Fmatrix=partialNumericDiff(states);

    Fmatrix=I+sampleTime*Fmatrix;
    
    states = states+sampleTime*predictEq';
    
    covP = (Fmatrix)*(covP*(Fmatrix')) + Q_noise;
   
    %updata=DDPH1,DDPH2,DDPH3,DDPH4,DDTH1,DDTH2,DDTH3,DDTH4,DPH1,DPH2,DPH3,DPH4,DSI1,DSI2,DTH1,DTH2,DTH3,DTH4,PH1,PH2,PH3,PH4,SI1,SI2,TH1,TH2,TH3,TH4
    updateEq=updateFunctionEKF(predictEq(3),predictEq(6),predictEq(8),predictEq(10),predictEq(2),predictEq(5),predictEq(7),predictEq(9),...
        states(3),states(6),states(8),states(10),states(1),states(4),states(2),states(5),states(7),states(9),...
        states(13),states(16),states(18),states(20),states(11),states(14),states(12),states(15),states(17),states(19));
    %H=DDPH1,DDPH2,DDPH3,DDPH4,DDTH1,DDTH2,DDTH3,DDTH4,DPH1,DPH2,DPH3,DPH4,DSI1,DSI2,DTH1,DTH2,DTH3,DTH4,PH1,PH2,PH3,PH4,SI1,SI2,TH1,TH2,TH3,TH4
    Hmatrix=HfunctionEKF(predictEq(3),predictEq(6),predictEq(8),predictEq(10),predictEq(2),predictEq(5),predictEq(7),predictEq(9),...
        states(3),states(6),states(8),states(10),states(1),states(4),states(2),states(5),states(7),states(9),...
        states(13),states(16),states(18),states(20),states(11),states(14),states(12),states(15),states(17),states(19));
    
    K=(covP*Hmatrix')/(Hmatrix*covP*Hmatrix'+R);
   
    z=[a1.signals.values(i,1) a1.signals.values(i,3) g1.signals.values(i,1) g1.signals.values(i,2) g1.signals.values(i,3) ma1.signals.values(i,1) ma1.signals.values(i,2) ma1.signals.values(i,3)...
       a2.signals.values(i,1) a2.signals.values(i,3) g2.signals.values(i,1) g2.signals.values(i,2) g2.signals.values(i,3) ma2.signals.values(i,1) ma2.signals.values(i,2) ma2.signals.values(i,3)...
       a3.signals.values(i,1) a3.signals.values(i,3) g3.signals.values(i,2) g3.signals.values(i,3) ma3.signals.values(i,1) ma3.signals.values(i,2) ma3.signals.values(i,3)...
       a4.signals.values(i,1) a4.signals.values(i,3) g4.signals.values(i,2) g4.signals.values(i,3) ma4.signals.values(i,1) ma4.signals.values(i,2) ma4.signals.values(i,3)]'; 

    states=states+K*(z-updateEq');
    covP=(I-K*Hmatrix)*covP;

    
    si1=states(11);
    th1=states(12);
    ph1=-states(13);
    si2=states(14);
    th2=states(15);
    ph2=-states(16);
    th3=states(17);
    ph3=-states(18);
    th4=states(19);
    ph4=-states(20);
    
    R_0_1=([[cos(ph1) -sin(ph1) 0];[sin(ph1) cos(ph1) 0];[0 0 1]]*[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]]*[[1 0 0];[0 cos(si1) -sin(si1)];[0 sin(si1) cos(si1)]]);
    R_0_2=([[cos(ph2) -sin(ph2) 0];[sin(ph2) cos(ph2) 0];[0 0 1]]*[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]]*[[1 0 0];[0 cos(si2) -sin(si2)];[0 sin(si2) cos(si2)]]);
    R_0_3=([[cos(ph3) -sin(ph3) 0];[sin(ph3) cos(ph3) 0];[0 0 1]]*[[cos(th3) 0 sin(th3)];[0 1 0];[-sin(th3) 0 cos(th3)]]);
    R_0_4=([[cos(ph4) -sin(ph4) 0];[sin(ph4) cos(ph4) 0];[0 0 1]]*[[cos(th4) 0 sin(th4)];[0 1 0];[-sin(th4) 0 cos(th4)]]);

    P1base=[0;0;0];
    P2tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0];
    P3tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[l3;0;0];
    P4tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[l3;0;0]+R_0_4*[l4;0;0];
    
    storedStates(i,:)=states(:)';
    storedPositions(i,:)=[P1base',P2tip',P3tip',P4tip'];
end

%%%%%%%%%%%%%%%cam data stuff 
costFunction=sqrt((sum((storedPositions(startTime:end,1)-camData.signals.values(startTime:end,1)).^2)+...
                   sum((storedPositions(startTime:end,2)-camData.signals.values(startTime:end,2)).^2)+...
                   sum((storedPositions(startTime:end,3)-camData.signals.values(startTime:end,3)).^2)+...
                   sum((storedPositions(startTime:end,4)-camData.signals.values(startTime:end,4)).^2)+...
                   sum((storedPositions(startTime:end,5)-camData.signals.values(startTime:end,5)).^2)+...
                   sum((storedPositions(startTime:end,6)-camData.signals.values(startTime:end,6)).^2)+...
                   sum((storedPositions(startTime:end,7)-camData.signals.values(startTime:end,7)).^2)+...
                   sum((storedPositions(startTime:end,8)-camData.signals.values(startTime:end,8)).^2)+...
                   sum((storedPositions(startTime:end,9)-camData.signals.values(startTime:end,9)).^2)+...
                   sum((storedPositions(startTime:end,10)-camData.signals.values(startTime:end,10)).^2)+...
                   sum((storedPositions(startTime:end,11)-camData.signals.values(startTime:end,11)).^2)+...
                   sum((storedPositions(startTime:end,12)-camData.signals.values(startTime:end,12)).^2))/(length(a1.signals.values(:,1))-startTime));
               
end

%
%
%
%