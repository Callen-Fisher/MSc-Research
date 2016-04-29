clc
clear
run('PoseCalc/runMySystem')
clc
syms th1 dth1 ddth1 th2 dth2 ddth2
%m1 m2 J1 J2 g l1 l2 J1 J2  A1 A2 Cd1 Cd2 rho magVecX magVecY magVecZ sampleTime Jr
disp('starting');
disp('calculating positions');


%the gen coordinates
q=[th1; th2];
dq=[dth1; dth2];
ddq=[ddth1; ddth2];

%the rotation matrices
Rpitch1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];

Rpitch2=[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];

R_0_1=simple(Rpitch1,'IgnoreAnalyticConstraints',true);
R_0_2=simple(Rpitch2,'IgnoreAnalyticConstraints',true);

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
T2=(sum(1/2*J1*dth1^2));
T2=T2+(sum(1/2*J2*dth2^2));
U=0;

disp('calculating gen forces');

Q=[0;0];

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
Dwpt=simple(systemEquations(2));
clear systemEquations M C G Q T1 T2 U Ttot Vtot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(' ');
disp('the Kalman Filter code');
disp(' ');
disp(' ');

disp('the F matrix');
states=[dth1 dth2 th1 th2];
Fequations=[Dwpm,Dwpt,dth1,dth2];
F = jacobian(Fequations,states);

disp('the measurement equations')
accMid=R_0_1.'*[0;0;g]+R_0_1.'*ddP1tip;
axMid=accMid(1);
ayMid=accMid(2);
azMid=accMid(3);
gyroMid=[0;dth1;0];
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
gyroTip=[0;dth2;0];
gyroXTip=gyroTip(1);
gyroYTip=gyroTip(2);
gyroZTip=gyroTip(3);
magTip=R_0_2.'*[magVecX;magVecY;magVecZ];
magXTip=magTip(1);
magYTip=magTip(2);
magZTip=magTip(3);

disp('the H matrix');
Hequations=[axMid,azMid,gyroYMid,magXMid,magZMid,axTip,azTip,gyroYTip,magXTip,magZTip];
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