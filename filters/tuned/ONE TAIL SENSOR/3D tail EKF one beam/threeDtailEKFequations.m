clc
clear
run('PoseCalc/runMySystem')
clc
syms th dth ddth ph dph  ddph 
disp('starting');
disp('calculating positions');


%the gen coordinates
q=[th; ph];
dq=[dth; dph];
ddq=[ddth; ddph];

%the inertia matrix
Jmid=[[0 0 0];[0 J1+J2 0];[0 0 J1+J2]];


%the rotation matrices
Rpitch=[[cos(th) 0 sin(th)];[0 1 0];[-sin(th) 0 cos(th)]];
Ryaw=[[cos(ph) -sin(ph) 0];[sin(ph) cos(ph) 0];[0 0 1]];

R=simple(Ryaw*Rpitch,'IgnoreAnalyticConstraints',true);


R_gyro=[[0 0 0];[0 cos(ph) 0];[0 0 1]];

%the positions
P1=R*[(l1+l2)/2; 0; 0];

P1tip=R*[l1+l2; 0; 0];

%the velocities
dP1=jacobian(P1,q)*dq;

dP1tip=jacobian(P1tip,q)*dq;

%the accelerations
ddP1tip=jacobian(dP1tip,q)*dq+jacobian(dP1tip,dq)*ddq;

disp('calculating the energies');
T1=(sum(1/2*(m1*transpose(dP1)*dP1)));
T2=(sum(1/2*Jmid*(R_gyro*[0;dth;dph]).^2));
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
Dwym=simple(systemEquations(2));
clear systemEquations M C G Q T1 T2 U Ttot Vtot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(' ');
disp('the Kalman Filter code');
disp(' ');
disp(' ');

disp('the F matrix');
states=[dth dph th ph];
Fequations=[Dwpm,Dwym,dth,dph];
F = jacobian(Fequations,states);

disp('the measurement equations')
accMid=R.'*[0;0;g]+R.'*ddP1tip;
axMid=accMid(1);
ayMid=accMid(2);
azMid=accMid(3);
gyroMid=R_gyro*[0;dth;dph];
gyroXMid=gyroMid(1);
gyroYMid=gyroMid(2);
gyroZMid=gyroMid(3);
magMid=R.'*[magVecX;magVecY;magVecZ];
magXMid=magMid(1);
magYMid=magMid(2);
magZMid=magMid(3);

disp('the H matrix');
Hequations=[axMid,azMid,gyroYMid,gyroZMid,magXMid,magYMid,magZMid];
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