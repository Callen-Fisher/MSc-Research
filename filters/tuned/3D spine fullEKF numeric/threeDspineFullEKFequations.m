clc
clear
run('PoseCalc/runMySystem')
clc
syms th1 dth1 ddth1 th2 dth2 ddth2 th3 dth3 ddth3 th4 dth4 ddth4
syms si1 dsi1 ddsi1 si2 dsi2 ddsi2 placeHolder signdph3 signdph4 signdth3 signdth4
syms ph1 dph1 ddph1 ph2 dph2 ddph2 ph3 dph3 ddph3 ph4 dph4 ddph4 Tp1 Tp2 Ty1 Ty2
%m1 m2 J1 J2 g l1 l2 J1 J2  A1 A2 Cd1 Cd2 rho magVecX magVecY magVecZ sampleTime Jr
disp('starting');
disp('calculating positions');

J1matrix=[[J1 0 0];[0 J1 0];[0 0 J1]];
J2matrix=[[J2 0 0];[0 J2 0];[0 0 J2]];
J3matrix=[[0 0 0];[0 J3 0];[0 0 J3]];
J4matrix=[[0 0 0];[0 J4 0];[0 0 J4]];
%the gen coordinates
q=[si1; th1; ph1; si2; th2; ph2; th3; ph3; th4; ph4];
dq=[dsi1; dth1; dph1; dsi2; dth2; dph2; dth3; dph3; dth4; dph4];
ddq=[ddsi1; ddth1; ddph1; ddsi2; ddth2; ddph2; ddth3; ddph3; ddth4; ddph4];

%the rotation matrices
Rroll1=[[1 0 0];[0 cos(si1) -sin(si1)];[0 sin(si1) cos(si1)]];
Rpitch1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
Ryaw1=[[cos(ph1) -sin(ph1) 0];[sin(ph1) cos(ph1) 0];[0 0 1]];

Rroll2=[[1 0 0];[0 cos(si2) -sin(si2)];[0 sin(si2) cos(si2)]];
Rpitch2=[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];
Ryaw2=[[cos(ph2) -sin(ph2) 0];[sin(ph2) cos(ph2) 0];[0 0 1]];

Rpitch3=[[cos(th3) 0 sin(th3)];[0 1 0];[-sin(th3) 0 cos(th3)]];
Ryaw3=[[cos(ph3) -sin(ph3) 0];[sin(ph3) cos(ph3) 0];[0 0 1]];

Rpitch4=[[cos(th4) 0 sin(th4)];[0 1 0];[-sin(th4) 0 cos(th4)]];
Ryaw4=[[cos(ph4) -sin(ph4) 0];[sin(ph4) cos(ph4) 0];[0 0 1]];

R_0_1=simple(Ryaw1*Rpitch1*Rroll1,'IgnoreAnalyticConstraints',true);
R_0_2=simple(Ryaw2*Rpitch2*Rroll2,'IgnoreAnalyticConstraints',true);
R_0_3=simple(Ryaw3*Rpitch3,'IgnoreAnalyticConstraints',true);
R_0_4=simple(Ryaw4*Rpitch4,'IgnoreAnalyticConstraints',true);
%the gyro rotation matrices
Rgyro_1_0=[[cos(ph1)*cos(th1) -sin(ph1) 0];[sin(ph1)*cos(th1) cos(ph1) 0];[-sin(th1) 0 1]];
Rgyro_2_0=[[cos(ph2)*cos(th2) -sin(ph2) 0];[sin(ph2)*cos(th2) cos(ph2) 0];[-sin(th2) 0 1]];
Rgyro_3_0=[[0 0 0];[0 cos(ph3) 0];[0 0 1]];
Rgyro_4_0=[[0 0 0];[0 cos(ph4) 0];[0 0 1]];
%the positions
P1=R_0_1*[l1/2; 0; 0];
P2=R_0_1*[l1; 0; 0]+R_0_2*[l2/2; 0; 0];
P3=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[l3/2;0;0];
P4=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[l3;0;0]+R_0_4*[l4/2;0;0];

P1base=[0;0;0];
P1tip=R_0_1*[l1; 0; 0];
P2tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0];
P3tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[l3;0;0];
P4tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[l3;0;0]+R_0_4*[l4;0;0];

%the velocities
dP1base=jacobian(P1base,q)*dq;
dP1=jacobian(P1,q)*dq;
dP2=jacobian(P2,q)*dq;
dP3=jacobian(P3,q)*dq;
dP4=jacobian(P4,q)*dq;

dP1tip=jacobian(P1tip,q)*dq;
dP2tip=jacobian(P2tip,q)*dq;
dP3tip=jacobian(P3tip,q)*dq;
dP4tip=jacobian(P4tip,q)*dq;
%the accelerations
ddP1base=jacobian(dP1base,q)*dq+jacobian(dP1base,dq)*ddq;
ddP1tip=jacobian(dP1tip,q)*dq+jacobian(dP1tip,dq)*ddq;
ddP2tip=jacobian(dP2tip,q)*dq+jacobian(dP2tip,dq)*ddq;
ddP3tip=jacobian(dP3tip,q)*dq+jacobian(dP3tip,dq)*ddq;
ddP4tip=jacobian(dP4tip,q)*dq+jacobian(dP4tip,dq)*ddq;

disp('calculating the energies');
T1=(sum(1/2*(m1*transpose(dP1)*dP1+m2*transpose(dP2)*dP2+m3*transpose(dP3)*dP3+m4*transpose(dP4)*dP4)));
T2=(sum(1/2*J1matrix*(Rgyro_1_0*[dsi1;dth1;dph1]).^2));
T2=T2+(sum(1/2*J2matrix*(Rgyro_2_0*[dsi2;dth2;dph2]).^2));
T2=T2+(sum(1/2*J3matrix*(Rgyro_3_0*[0;dth3;dph3]).^2));
T2=T2+(sum(1/2*J4matrix*(Rgyro_4_0*[0;dth4;dph4]).^2));
U=m3*g*P3(3)+m4*g*P4(3);

disp('calculating gen forces');

Q=[placeHolder;placeHolder;placeHolder;placeHolder;placeHolder;placeHolder;placeHolder;placeHolder;placeHolder;placeHolder];

dy1=1/2*Cd1*rho*A1*dph3^2*l3^2*signdph3;
dp1=1/2*Cd1*rho*A1*(dth3*cos(ph3))^2*l3^2*signdth3;
dy2=1/2*Cd2*rho*A2*(dph4^2*signdph4)*l4^2;
dp2=1/2*Cd2*rho*A2*(dth4*cos(ph4))^2*signdth4*l4^2;
 
f1=R_0_3*[0;-dy1;dp1];
f2=R_0_4*[0;-dy2;dp2];
 
r1=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[3*l3/4;0;0];
r2=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0]+R_0_3*[l3;0;0]+R_0_4*[3*l4/4;0;0];

f=[f1,f2];
r=[r1,r2];
dr=symMat([length(r(:,1)) length(q)*length(r(1,:))],'dr','real');
counter=1;
for j=1:1:length(q)
    for i=1:1:length(f(1,:))
        if(Q(j)==placeHolder)
            dr(:,counter)=diff(r(:,i),q(j));
            Q(j)=sum(f(:,i).*dr(:,counter));
            counter=counter+1;
        else
            dr(:,counter)=diff(r(:,i),q(j));
            Q(j)=Q(j)+sum(f(:,i).*dr(:,counter));
            counter=counter+1;
        end
    end
end
Q=[0;0;0;0;0;0;Q(7)+Tp1;Q(8)+Ty1;Q(9)+Tp2;Q(10)+Ty2];


disp('simplify the energies');
Ttot = T1+T2;
Ttot = simple(Ttot,'IgnoreAnalyticConstraints',true);
Vtot = U;
% Mass Matrix
disp('calculating the M matrix')
M=jacobian(jacobian(Ttot,dq).',dq);
M = simple(M,'IgnoreAnalyticConstraints',true);
% C Matrix -- contains the centrifugal and coriolis accelerations
disp('calculating the C matrix')
C = symMat([length(q) length(q)],'C','real');
for i = 1 : length(q)
    for j = 1 : length(q)
        for k = 1: length(q)
            if k == 1
                C(i,j) = 0.5*(diff(M(i,j),q(k)) + diff(M(i,k),q(j)) - diff(M(j,k),q(i)))*dq(k);
            else
                C(i,j) = C(i,j)+ 0.5*(diff(M(i,j),q(k)) + diff(M(i,k),q(j)) - diff(M(j,k),q(i)))*dq(k);
            end
        end        
    end
end
C = simple(C,'IgnoreAnalyticConstraints',true);
% G Matrix --> Contains the potential energy
disp('calculating the G matrix')
G = symMat([length(q) 1],'g','real');
for i = 1: length(q)
    G(i) = (diff(Vtot,q(i)));
end
G = simple(G,'IgnoreAnalyticConstraints',true);

%disp('calculating the system equations')
%systemEquations=((M)\(-C*dq - G + Q));

% disp('solving equations');
% disp('1');
% Dwsi1=simple(systemEquations(1));
% disp('2');
% Dwth1=simple(systemEquations(2));
% disp('3');
% Dwph1=simple(systemEquations(3));
% disp('4');
% Dwsi2=simple(systemEquations(4));
% disp('5');
% Dwth2=simple(systemEquations(5));
% disp('6');
% Dwph2=simple(systemEquations(6));
% disp('7');
% Dwth3=simple(systemEquations(7));
% disp('8');
% Dwph3=simple(systemEquations(8));
% disp('9');
% Dwth4=simple(systemEquations(9));
% disp('10');
% Dwph4=simple(systemEquations(10));
% clear systemEquations M C G Q T1 T2 U Ttot Vtot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(' ');
disp('the Kalman Filter code');
disp(' ');
disp(' ');

disp('the F matrix');
states=[Tp1 Ty1 Tp2 Ty2 dsi1 dth1 dph1 dsi2 dth2 dph2 dth3 dph3 dth4 dph4 si1 th1 ph1 si2 th2 ph2 th3 ph3 th4 ph4];
%Fequations=[Dwsi1 Dwth1 Dwph1 Dwsi2 Dwth2 Dwph2 Dwth3 Dwph3 Dwth4 Dwph4 dsi1 dth1 dph1 dsi2 dth2 dph2 dth3 dph3 dth4 dph4];
%F = jacobian(Fequations,states);

disp('the measurement equations')
acc1=R_0_1.'*[0;0;g]+R_0_1.'*ddP1base;
ax1=acc1(1);
ay1=acc1(2);
az1=acc1(3);
gyro1=Rgyro_1_0*[dsi1;dth1;dph1];
gyroX1=gyro1(1);
gyroY1=gyro1(2);
gyroZ1=gyro1(3);
mag1=R_0_1.'*[magVecX;magVecY;magVecZ];
magX1=mag1(1);
magY1=mag1(2);
magZ1=mag1(3);

acc2=R_0_2.'*[0;0;g]+R_0_2.'*ddP2tip;
ax2=acc2(1);
ay2=acc2(2);
az2=acc2(3);
gyro2=Rgyro_2_0*[dsi2;dth2;dph2];
gyroX2=gyro2(1);
gyroY2=gyro2(2);
gyroZ2=gyro2(3);
mag2=R_0_2.'*[magVecX;magVecY;magVecZ];
magX2=mag2(1);
magY2=mag2(2);
magZ2=mag2(3);

acc3=R_0_3.'*[0;0;g]+R_0_3.'*ddP3tip;
ax3=acc3(1);
ay3=acc3(2);
az3=acc3(3);
gyro3=Rgyro_3_0*[0;dth3;dph3];
gyroX3=gyro3(1);
gyroY3=gyro3(2);
gyroZ3=gyro3(3);
mag3=R_0_3.'*[magVecX;magVecY;magVecZ];
magX3=mag3(1);
magY3=mag3(2);
magZ3=mag3(3);

acc4=R_0_4.'*[0;0;g]+R_0_4.'*ddP4tip;
ax4=acc4(1);
ay4=acc4(2);
az4=acc4(3);
gyro4=Rgyro_4_0*[0;dth4;dph4];
gyroX4=gyro4(1);
gyroY4=gyro4(2);
gyroZ4=gyro4(3);
mag4=R_0_4.'*[magVecX;magVecY;magVecZ];
magX4=mag4(1);
magY4=mag4(2);
magZ4=mag4(3);

disp('the H matrix');
Hequations=[ax1,az1,gyroX1,gyroY1,gyroZ1,magX1,magY1,magZ1...
           ,ax2,az2,gyroX2,gyroY2,gyroZ2,magX2,magY2,magZ2...
           ,ax3,az3,gyroY3,gyroZ3,magX3,magY3,magZ3...
           ,ax4,az4,gyroY4,gyroZ4,magX4,magY4,magZ4];
H = jacobian(Hequations,states);

%predict=Fequations;
update=Hequations;

%C,M,G,Q     
save('EKFmatrices/CFullEKF','C');
save('EKFmatrices/MFullEKF','M');
save('EKFmatrices/GFullEKF','G');
% save('EKFmatrices/QFullEKF','Q');
save('EKFmatrices/HFullEKF','H');
save('EKFmatrices/updateFullEKF','update');
save('EKFmatrices/statesFullEKF','states');


disp('matlab function C')
matlabFunction(C,'file','CfunctionFullEKF');
disp('matlab function M')
matlabFunction(M,'file','MfunctionFullEKF');
disp('matlab function G')
matlabFunction(G,'file','GfunctionFullEKF');
disp('matlab function Q')
matlabFunction(Q,'file','QfunctionFullEKF');
disp('matlab function H')
matlabFunction(H,'file','HfunctionFullEKF');
disp('matlab function update')
matlabFunction(update,'file','updateFunctionFullEKF');

clear all
clc