clc
clear
run('PoseCalc/runMySystem')
clc
syms th1 dth1 ddth1 th2 dth2 ddth2 th3 dth3 ddth3 th4 dth4 ddth4
%m1 m2 J1 J2 g l1 l2 J1 J2  A1 A2 Cd1 Cd2 rho magVecX magVecY magVecZ sampleTime Jr
disp('starting');
disp('calculating positions');


%the gen coordinates
q=[th1; th2; th3; th4];
dq=[dth1; dth2; dth3; dth4];
ddq=[ddth1; ddth2; ddth3; ddth4];

%the rotation matrices
Rpitch1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
Rpitch2=[[cos(th2) 0 sin(th2)];[0 1 0];[-sin(th2) 0 cos(th2)]];
Rpitch3=[[cos(th3) 0 sin(th3)];[0 1 0];[-sin(th3) 0 cos(th3)]];
Rpitch4=[[cos(th4) 0 sin(th4)];[0 1 0];[-sin(th4) 0 cos(th4)]];

R_0_1=simple(Rpitch1,'IgnoreAnalyticConstraints',true);
R_0_2=simple(Rpitch2,'IgnoreAnalyticConstraints',true);
R_0_3=simple(Rpitch3,'IgnoreAnalyticConstraints',true);
R_0_4=simple(Rpitch4,'IgnoreAnalyticConstraints',true);

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
T2=(sum(1/2*J1*dth1^2));
T2=T2+(sum(1/2*J2*dth2^2));
T2=T2+(sum(1/2*J3*dth3^2));
T2=T2+(sum(1/2*J4*dth4^2));
U=0;

disp('calculating gen forces');

Q=[0;0;0;0];

disp('simplify the energies');
Ttot = T1+T2;
Ttot = simple(Ttot,'IgnoreAnalyticConstraints',true);
Vtot = U;
% Mass Matrix
disp('calculating the M matrix')
M=jacobian(jacobian(Ttot,dq).',dq)
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
Dwp1=simple(systemEquations(1));
disp('2');
Dwp2=simple(systemEquations(2));
disp('3');
Dwp3=simple(systemEquations(3));
disp('4');
Dwp4=simple(systemEquations(4));
clear systemEquations M C G Q T1 T2 U Ttot Vtot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(' ');
disp('the Kalman Filter code');
disp(' ');
disp(' ');

disp('the F matrix');
states=[dth1 dth2 dth3 dth4 th1 th2 th3 th4];
Fequations=[Dwp1,Dwp2,Dwp3,Dwp4,dth1,dth2,dth3,dth4];
F = jacobian(Fequations,states);

disp('the measurement equations')
acc1=R_0_1.'*[0;0;g]+R_0_1.'*ddP1base;
ax1=acc1(1);
ay1=acc1(2);
az1=acc1(3);
gyro1=[0;dth1;0];
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
gyro2=[0;dth2;0];
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
gyro3=[0;dth3;0];
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
gyro4=[0;dth4;0];
gyroX4=gyro4(1);
gyroY4=gyro4(2);
gyroZ4=gyro4(3);
mag4=R_0_4.'*[magVecX;magVecY;magVecZ];
magX4=mag4(1);
magY4=mag4(2);
magZ4=mag4(3);

disp('the H matrix');
Hequations=[ax1,az1,gyroY1,magX1,magZ1,ax2,az2,gyroY2,magX2,magZ2,ax3,az3,gyroY3,magX3,magZ3,ax4,az4,gyroY4,magX4,magZ4];
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