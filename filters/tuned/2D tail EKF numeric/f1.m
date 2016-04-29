function [ val ] = f1( states,epsilon )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Q=[0;0];
M=MfunctionEKF(states(3),states(4));
C=CfunctionEKF(states(1)+epsilon,states(2),states(3),states(4));
dq=[states(1)+epsilon; states(2)];
G=[0;0];
systemEquations1=(M)\(-C*dq - G + Q);
Q=[0;0];
M=MfunctionEKF(states(3),states(4));
C=CfunctionEKF(states(1)-epsilon,states(2),states(3),states(4));
dq=[states(1)-epsilon; states(2)];
G=[0;0];
systemEquations2=(M)\(-C*dq - G + Q);
val1_1=(systemEquations1(1)-systemEquations2(1))/(2*epsilon);
val1_2=(systemEquations1(2)-systemEquations2(2))/(2*epsilon);

Q=[0;0];
M=MfunctionEKF(states(3),states(4));
C=CfunctionEKF(states(1),states(2)+epsilon,states(3),states(4));
dq=[states(1); states(2)+epsilon];
G=[0;0];
systemEquations1=(M)\(-C*dq - G + Q);
Q=[0;0];
M=MfunctionEKF(states(3),states(4));
C=CfunctionEKF(states(1),states(2)-epsilon,states(3),states(4));
dq=[states(1); states(2)-epsilon];
G=[0;0];
systemEquations2=(M)\(-C*dq - G + Q);
val2_1=(systemEquations1(1)-systemEquations2(1))/(2*epsilon);
val2_2=(systemEquations1(2)-systemEquations2(2))/(2*epsilon);

Q=[0;0];
M=MfunctionEKF(states(3)+epsilon,states(4));
C=CfunctionEKF(states(1),states(2),states(3)+epsilon,states(4));
dq=[states(1); states(2)];
G=[0;0];
systemEquations1=(M)\(-C*dq - G + Q);
Q=[0;0];
M=MfunctionEKF(states(3)-epsilon,states(4));
C=CfunctionEKF(states(1),states(2),states(3)-epsilon,states(4));
dq=[states(1); states(2)];
G=[0;0];
systemEquations2=(M)\(-C*dq - G + Q);
val3_1=(systemEquations1(1)-systemEquations2(1))/(2*epsilon);
val3_2=(systemEquations1(2)-systemEquations2(2))/(2*epsilon);

Q=[0;0];
M=MfunctionEKF(states(3),states(4)+epsilon);
C=CfunctionEKF(states(1),states(2),states(3),states(4)+epsilon);
dq=[states(1); states(2)];
G=[0;0];
systemEquations1=(M)\(-C*dq - G + Q);
Q=[0;0];
M=MfunctionEKF(states(3),states(4)-epsilon);
C=CfunctionEKF(states(1),states(2),states(3),states(4)-epsilon);
dq=[states(1); states(2)];
G=[0;0];
systemEquations2=(M)\(-C*dq - G + Q);
val4_1=(systemEquations1(1)-systemEquations2(1))/(2*epsilon);
val4_2=(systemEquations1(2)-systemEquations2(2))/(2*epsilon);

val=[[val1_1,val2_1,val3_1,val4_1];[val1_2,val2_2,val3_2,val4_2]];





end

