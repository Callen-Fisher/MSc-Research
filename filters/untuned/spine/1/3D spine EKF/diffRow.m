function [ val ] = diffRow( states1,states2,e )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    Q=[0;0;0;0;0;0;0;0;0;0];
    M=MfunctionEKF(states1(13),states1(16),states1(18),states1(20),states1(12),states1(15),states1(17),states1(19));
    C=CfunctionEKF(states1(3),states1(6),states1(8),states1(10),states1(1),states1(4),states1(2),states1(5),states1(7),states1(9),...
        states1(13),states1(16),states1(18),states1(20),states1(12),states1(15),states1(17),states1(19));
    dq=[states1(1); states1(2); states1(3); states1(4); states1(5); states1(6); states1(7); states1(8); states1(9); states1(10)];
    G=[0;0;0;0;0;0;0;0;0;0];
    systemEquations1=(M)\(-C*dq - G + Q);
    Q=[0;0;0;0;0;0;0;0;0;0];
    M=MfunctionEKF(states2(13),states2(16),states2(18),states2(20),states2(12),states2(15),states2(17),states2(19));
    C=CfunctionEKF(states2(3),states2(6),states2(8),states2(10),states2(1),states2(4),states2(2),states2(5),states2(7),states2(9),...
        states2(13),states2(16),states2(18),states2(20),states2(12),states2(15),states2(17),states2(19));
    dq=[states2(1); states2(2); states2(3); states2(4); states2(5); states2(6); states2(7); states2(8); states2(9); states2(10)];
    G=[0;0;0;0;0;0;0;0;0;0];
    systemEquations2=(M)\(-C*dq - G + Q);
    
    
    
    val=(systemEquations1-systemEquations2)/(2*e);
end

