clc
clear
run('PoseCalc/runMySystem')
clc
syms th1 dth1 ddth1 th2 dth2 ddth2 ph1 dph1  ddph1 ph2 dph2 ddph2 
%m1 m2 J1 J2 g l1 l2 J1 J2  A1 A2 Cd1 Cd2 rho magVecX magVecY magVecZ sampleTime Jr
disp('starting');
disp('calculating positions');


%the gen coordinates
q=[th1; ph1; th2; ph2];
dq=[dth1; dph1; dth2; dph2];
ddq=[ddth1; ddph1; ddth2; ddph2];

%the inertia matrix
Jmid=[[0 0 0];[0 J1 0];[0 0 J1]];
Jtip=[[0 0 0];[0 J2 0];[0 0 J2]];

%the rotation matrices
Rpitch1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
Ryaw1=[[cos(ph1) -sin(ph1) 0];[sin(ph1) cos(ph1) 0];[0 0 1]];

Rpitch2=[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];
Ryaw2=[[cos(ph2) -sin(ph2) 0];[sin(ph2) cos(ph2) 0];[0 0 1]];

R_0_1=simple(Ryaw1*Rpitch1,'IgnoreAnalyticConstraints',true);
R_0_2=simple(Ryaw2*Rpitch2,'IgnoreAnalyticConstraints',true);

R_gyro_1_0=[[0 0 0];[0 cos(ph1) 0];[0 0 1]];
R_gyro_2_1=[[0 0 0];[0 cos(ph2) 0];[0 0 1]];

%the positions
P1=R_0_1*[l1/2; 0; 0];
P2=R_0_1*[l1; 0; 0]+R_0_2*[l2/2; 0; 0];

P1tip=R_0_1*[l1; 0; 0];
P2tip=R_0_1*[l1; 0; 0]+R_0_2*[l2; 0; 0];

%the velocities
dP1=jacobian(P1,q)*dq;
dP2=jacobian(P2,q)*dq;

dP1tip=jacobian(P1tip,q)*dq;
dP2tip=jacobian(P2tip,q)*dq;

%the accelerations
ddP1tip=jacobian(dP1tip,q)*dq+jacobian(dP1tip,dq)*ddq;
ddP2tip=jacobian(dP2tip,q)*dq+jacobian(dP2tip,dq)*ddq;

disp('calculating the energies');
T1=(sum(1/2*(m1*transpose(dP1)*dP1+m2*transpose(dP2)*dP2)));
T2=(sum(1/2*Jmid*(R_gyro_1_0*[0;dth1;dph1]).^2));
T2=T2+(sum(1/2*Jtip*(R_gyro_2_1*[0;dth2;dph2]).^2));
U=0;

disp('calculating gen forces');

Q=[0;0;0;0];

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

disp('calculating the system equations')
systemEquations=((M)\(-C*dq - G + Q));

disp('solving equations');
disp('1');
Dwpm=simple(systemEquations(1));
disp('2');
Dwym=simple(systemEquations(2));
disp('3');
Dwpt=simple(systemEquations(3));
disp('4');
Dwyt=simple(systemEquations(4));
clear systemEquations M C G Q T1 T2 U Ttot Vtot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(' ');
disp('the Kalman Filter code');
disp(' ');
disp(' ');

disp('the F matrix');
states=[dth1 dph1 dth2 dph2 th1 ph1 th2 ph2];
Fequations=[Dwpm,Dwym,Dwpt,Dwyt,dth1,dph1,dth2,dph2];
F = jacobian(Fequations,states);

disp('the measurement equations')
accMid=R_0_1.'*[0;0;g]+R_0_1.'*ddP1tip;
axMid=accMid(1);
ayMid=accMid(2);
azMid=accMid(3);
gyroMid=R_gyro_1_0*[0;dth1;dph1];
gyroXMid=gyroMid(1);
gyroYMid=gyroMid(2);
gyroZMid=gyroMid(3);
magMid=R_0_1.'*[magVecX;magVecY;magVecZ];
magXMid=magMid(1);
magYMid=magMid(2);
magZMid=magMid(3);

accTip=R_0_2.'*[0;0;g]+R_0_2.'*ddP2tip;
axTip=accTip(1);
ayTip=accTip(2);
azTip=accTip(3);
gyroTip=R_gyro_2_1*[0;dth2;dph2];
gyroXTip=gyroTip(1);
gyroYTip=gyroTip(2);
gyroZTip=gyroTip(3);
magTip=R_0_2.'*[magVecX;magVecY;magVecZ];
magXTip=magTip(1);
magYTip=magTip(2);
magZTip=magTip(3);

disp('the H matrix');
Hequations=[axMid,azMid,gyroYMid,gyroZMid,magXMid,magYMid,magZMid,axTip,azTip,gyroYTip,gyroZTip,magXTip,magYTip,magZTip];
H = jacobian(Hequations,states);

predict=Fequations;
update=Hequations;


save('EKFmatrices/HEKF','H');
save('EKFmatrices/FEKF','F');
save('EKFmatrices/predictEKF','predict');
save('EKFmatrices/updateEKF','update');
save('EKFmatrices/statesEKF','states');

disp('matlab function F')
matlabFunction(F,'file','FfunctionEKF');
disp('matlab function H')
matlabFunction(H,'file','HfunctionEKF');
disp('matlab function predict')
matlabFunction(predict,'file','predictFunctionEKF');
disp('matlab function update')
matlabFunction(update,'file','updateFunctionEKF');

clear all
clc