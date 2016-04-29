clc
clear
run('PoseCalc/runMySystem')
clc
syms th dth ddth ph dph ddph
syms Tp Ty
syms signdth signdph placeHolder
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
U=(m1+m2)*g*P1(3);

disp('calculating gen forces');

Q=[placeHolder;placeHolder];

dy1=1/2*Cd1*rho*(A1+A2)*dph^2*(l1+l2)^2*signdph;
dp1=1/2*Cd1*rho*(A1+A2)*(dth*cos(ph))^2*(l1+l2)^2*signdth;
 
f1=R*[0;-dy1;dp1];
 
r1=R*[3*(l1+l2)/4;0;0];

f=[f1];
r=[r1];
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
Q=[Q(1)+Tp;Q(2)+Ty];

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
Dwym=simple(systemEquations(2));
clear systemEquations M C G Q T1 T2 U Ttot Vtot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(' ');
disp('the Kalman Filter code');
disp(' ');
disp(' ');

DTpm=0;
DTym=0;
DTpt=0;
DTyt=0;

disp('the F matrix');
states=[Tp Ty dth dph th ph];
Fequations=[DTpm,DTym,Dwpm,Dwym,dth,dph];
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