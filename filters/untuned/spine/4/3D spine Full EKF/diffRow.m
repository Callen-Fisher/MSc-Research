function [ val ] = diffRow( states1,states2,e )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



    %Q=TP1,TP2,TY1,TY2,DPH3,DPH4,DTH3,DTH4,PH3,PH4,SIGNDPH3,SIGNDPH4,SIGNDTH3,SIGNDTH4,TH3,TH4
    Q=QfunctionFullEKF(states1(1),states1(3),states1(2),states1(4),states1(12),states1(14),states1(11),states1(13),states1(22),states1(24),...
        sign(states1(12)),sign(states1(14)),sign(states1(11)),sign(states1(13)),states1(21),states1(23));%was an error
    %M=PH1,PH2,PH3,PH4,TH1,TH2,TH3,TH4
    M=MfunctionFullEKF(states1(17),states1(20),states1(22),states1(24),states1(16),states1(19),states1(21),states1(23));%ERROR
    %C=DPH1,DPH2,DPH3,DPH4,DSI1,DSI2,DTH1,DTH2,DTH3,DTH4,PH1,PH2,PH3,PH4,TH1,TH2,TH3,TH4
    C=CfunctionFullEKF(states1(7),states1(10),states1(12),states1(14),states1(5),states1(8),states1(6),states1(9),states1(11),states1(13),states1(17),states1(20),states1(22),states1(24),states1(16),states1(19),states1(21),states1(23));%ERROR
    dq=[states1(5); states1(6); states1(7); states1(8); states1(9); states1(10); states1(11); states1(12); states1(13); states1(14)];
    %G=TH1,TH2,TH3,TH4
    G=GfunctionFullEKF(states1(16),states1(19),states1(21),states1(23));
    systemEquations1=(M)\(-C*dq - G + Q);

    %Q=TP1,TP2,TY1,TY2,DPH3,DPH4,DTH3,DTH4,PH3,PH4,SIGNDPH3,SIGNDPH4,SIGNDTH3,SIGNDTH4,TH3,TH4
    Q=QfunctionFullEKF(states2(1),states2(3),states2(2),states2(4),states2(12),states2(14),states2(11),states2(13),states2(22),states2(24),...
        sign(states2(12)),sign(states2(14)),sign(states2(11)),sign(states2(13)),states2(21),states2(23));%was an error
    %M=PH1,PH2,PH3,PH4,TH1,TH2,TH3,TH4
    M=MfunctionFullEKF(states2(17),states2(20),states2(22),states2(24),states2(16),states2(19),states2(21),states2(23));%ERROR
    %C=DPH1,DPH2,DPH3,DPH4,DSI1,DSI2,DTH1,DTH2,DTH3,DTH4,PH1,PH2,PH3,PH4,TH1,TH2,TH3,TH4
    C=CfunctionFullEKF(states2(7),states2(10),states2(12),states2(14),states2(5),states2(8),states2(6),states2(9),states2(11),states2(13),states2(17),states2(20),states2(22),states2(24),states2(16),states2(19),states2(21),states2(23));%ERROR
    dq=[states2(5); states2(6); states2(7); states2(8); states2(9); states2(10); states2(11); states2(12); states2(13); states2(14)];
    %G=TH1,TH2,TH3,TH4
    G=GfunctionFullEKF(states2(16),states2(19),states2(21),states2(23));
    systemEquations2=(M)\(-C*dq - G + Q);


    val=(systemEquations1-systemEquations2)/(2*e);
  
end

