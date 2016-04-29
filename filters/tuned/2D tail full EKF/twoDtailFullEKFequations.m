clc
clear
run('PoseCalc/runMySystem')
clc
syms th1 dth1 ddth1 th2 dth2 ddth2 
syms Tp1 Tp2 
%m1 m2 J1 J2 g l1 l2 J1 J2  A1 A2 Cd1 Cd2 rho magVecX magVecY magVecZ sampleTime Jr
syms signdth1 signdth2 placeHolder
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
U=m1*g*P1(3)+m2*g*P2(3);

disp('calculating gen forces');

Q=[placeHolder;placeHolder];

dp1=1/2*Cd1*rho*A1*(dth1)^2*l1^2*signdth1;
dp2=1/2*Cd2*rho*A2*(dth2)^2*l2^2*signdth2;
 
f1=R_0_1*[0;0;dp1];
f2=R_0_2*[0;0;dp2];
 
r1=R_0_1*[3*l1/4;0;0];
r2=R_0_1*[l1;0;0]+R_0_2*[3*l2/4;0;0];

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
Q=[Q(1)+Tp1;Q(2)+Tp2];

disp('simplify the energies');
Ttot = T1+T2;
Ttot = simple(Ttot,'IgnoreAnalyticConstraints',true);
Vtot = U;
Vtot = simple(Vtot,'IgnoreAnalyticConstraints',true);
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

DTpm=0;
DTpt=0;

disp('the F matrix');
states=[Tp1 Tp2 dth1 dth2 th1 th2];
Fequations=[DTpm,DTpt,Dwpm,Dwpt,dth1,dth2];
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


save('fullEKFmatrices/HfullEKF','H');
save('fullEKFmatrices/FfullEKF','F');
save('fullEKFmatrices/predictFullEKF','predict');
save('fullEKFmatrices/updateFullEKF','update');
save('fullEKFmatrices/statesFullEKF','states');

disp('matlab function F')
matlabFunction(F,'file','FfunctionFullEKF');
disp('matlab function H')
matlabFunction(H,'file','HfunctionFullEKF');
disp('matlab function predict')
matlabFunction(predict,'file','predictFunctionFullEKF');
disp('matlab function update')
matlabFunction(update,'file','updateFunctionFullEKF');

clear all
clc