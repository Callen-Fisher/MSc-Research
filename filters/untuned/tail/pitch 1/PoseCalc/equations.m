clc
clear

syms ph1 th1 ph2 th2 l1 l2 
disp('starting');
disp('calculating positions');

R1wrt0=[[cos(ph1)*cos(th1) -sin(ph1) cos(ph1)*sin(th1)];...
        [sin(ph1)*cos(th1) cos(ph1)  sin(ph1)*sin(th1)];...
        [-sin(th1)         0         cos(th1)]];%the rotation matrix
R2wrt1=[[cos(ph2)*cos(th2) -sin(ph2) cos(ph2)*sin(th2)];...
        [sin(ph2)*cos(th2) cos(ph2)  sin(ph2)*sin(th2)];...
        [-sin(th2)         0         cos(th2)]];
R2wrt0=R1wrt0*R2wrt1;

R1wrt0T=[[cos(ph1)*cos(th1) sin(ph1)*cos(th1)   -sin(th1)];...
        [-sin(ph1)          cos(ph1)            0];...
        [cos(ph1)*sin(th1)  sin(ph1)*sin(th1)   cos(th1)]];%the rotation matrix
R2wrt1T=[[cos(ph2)*cos(th2)  sin(ph2)*cos(th2)  -sin(th2)];...
        [-sin(ph2)           cos(ph2)           0];...
        [cos(ph2)*sin(th2)   sin(ph2)*sin(th2)  cos(th2)]];
R2wrt0T=R1wrt0T*R2wrt1T;

Rgyro1wrt0=[[1 0];[0 cos(th1)]];%the rotation matrix for the gyro
Rgyro2wrt1=[[1 0];[0 cos(th2)]];
Rgyro2wrt0=Rgyro1wrt0*Rgyro2wrt1;

Rgyro1wrt0T=[[1 0];[0 cos(th1)]];%the rotation matrix for the gyro
Rgyro2wrt1T=[[1 0];[0 cos(th2)]];
Rgyro2wrt0T=Rgyro1wrt0T*Rgyro2wrt1T;

P1=R1wrt0*[l1/2; 0; 0];%position of mass 1 (in the middle of beam 1)
P2=R1wrt0*[l1; 0; 0]+R2wrt0*[l2/2; 0; 0];%position of mass 2 (in the middle of beam 2)

dP1=simple(diff(P1,'th1')*dth1+diff(P1,'th2')*dth2+diff(P1,'ph1')*dph1+diff(P1,'ph2')*dph2);%the velocities of the masses
dP2=simple(diff(P2,'th1')*dth1+diff(P2,'th2')*dth2+diff(P2,'ph1')*dph1+diff(P2,'ph2')*dph2);%wrt time so multiply by d..  (dth1 etc)
disp('calculating energies');
T1=sum(1/2*(m1*dP1.^2+m2*dP2.^2));%kinetic energy due to translation
T2=1/2*(J1*(dth1^2+dph1^2)+J2*(dth2^2+dph2^2));%kinetic energy due to rotation 
U=m1*g*P1(3)+m2*g*P2(3);%potential energy of the masses

L=T1+T2-U;%the lagrange equation 
disp('differentiating 1 of 3');
dL_dx=simple([diff(L,'th1');diff(L,'ph1');diff(L,'th2');diff(L,'ph2')]);
disp('differentiating 2 of 3');
dL_ddx=simple([diff(L,'dth1');diff(L,'dph1');diff(L,'dth2');diff(L,'dph2')]);
disp('differentiating 3 of 3');
d_dt_dL_ddx=diff(dL_ddx,'th1')*dth1+diff(dL_ddx,'dth1')*ddth1+diff(dL_ddx,'th2')*dth2+diff(dL_ddx,'dth2')*ddth2+diff(dL_ddx,'ph1')*dph1+diff(dL_ddx,'dph1')*ddph1+diff(dL_ddx,'ph2')*dph2+diff(dL_ddx,'dph2')*ddph2;

disp('calculating gen forces');
F1=-1/2*rho*Cd1*A1*l1^2*dth1^2;%the forces due to drag
F2=-1/2*rho*Cd1*A1*l1^2*dph1^2;
F3=-1/2*rho*Cd2*A2*l2^2*dth2^2;
F4=-1/2*rho*Cd2*A2*l2^2*dph2^2;
F5=-1/2*rho*Cd2*A2*cos(th2)^3*l2^2*dth1^2;%forces due to second beam rotating by th1 and ph1
F6=-1/2*rho*Cd2*A2*cos(ph2)^3*l2^2*dph1^2;

R1=3*l1/4;%the radius of where the forces act 
R2=3*l1/4;
R3=3*l2/4;
R4=3*l2/4;
R5=l1+3*l2*cos(th2)/4;
R6=l1+3*l2*cos(ph2)/4;

disp('calculating drag 1 of 4');
dragQth1=F1*R1+F5*R5;%the drag torques
disp('calculating drag 2 of 4');
dragQph1=F2*R2+F6*R6;
disp('calculating drag 3 of 4');
dragQth2=F3*R3;
disp('calculating drag 4 of 4');
dragQph2=F4*R4;

Q=[Tp1+dragQth1;Ty1+dragQph1;Tp2+dragQth2;Ty2+dragQph2];%the generalized torques
disp('calculating equations');
systemEquations=simple(dL_dx-d_dt_dL_ddx-Q);
disp('solving equations');
ddth1Eq=simple(solve(systemEquations(1),'ddth1'))%the equations 
ddph1Eq=simple(solve(systemEquations(2),'ddph1'))
ddth2Eq=simple(solve(systemEquations(3),'ddth2'))
ddph2Eq=simple(solve(systemEquations(4),'ddph2'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp(' ');
disp('the Kalman Filter code');
disp(' ');
disp(' ');



disp('the F matrix');
disp('calculating F 1 of 12');
F1=[0                          0                         0                          0                         0                           0                           0                           0                           0                          0                          0                         0                          ];
disp('calculating F 2 of 12');
F2=[0                          0                         0                          0                         0                           0                           0                           0                           0                          0                          0                         0                          ];
disp('calculating F 3 of 12');
F3=[0                          0                         0                          0                         0                           0                           0                           0                           0                          0                          0                         0                          ];
disp('calculating F 4 of 12');
F4=[0                          0                         0                          0                         0                           0                           0                           0                           0                          0                          0                         0                          ];
disp('calculating F 5 of 12');
F5=[simple(diff(Dwpm,'Tp1'))   simple(diff(Dwpm,'Ty1'))  simple(diff(Dwpm,'Tp2'))   simple(diff(Dwpm,'Ty2'))  simple(diff(Dwpm,'dth1'))   simple(diff(Dwpm,'dph1'))  simple(diff(Dwpm,'dth2'))   simple(diff(Dwpm,'dph2'))  simple(diff(Dwpm,'th1'))   simple(diff(Dwpm,'ph1'))  simple(diff(Dwpm,'th2'))  simple(diff(Dwpm,'ph2'))  ];
disp('calculating F 6 of 12');
F6=[simple(diff(Dwym,'Tp1'))   simple(diff(Dwym,'Ty1'))  simple(diff(Dwym,'Tp2'))   simple(diff(Dwym,'Ty2'))  simple(diff(Dwym,'dth1'))   simple(diff(Dwym,'dph1'))  simple(diff(Dwym,'dth2'))   simple(diff(Dwym,'dph2'))  simple(diff(Dwym,'th1'))   simple(diff(Dwym,'ph1'))  simple(diff(Dwym,'th2'))  simple(diff(Dwym,'ph2'))  ];
disp('calculating F 7 of 12');
F7=[simple(diff(Dwpt,'Tp1'))   simple(diff(Dwpt,'Ty1'))  simple(diff(Dwpt,'Tp2'))   simple(diff(Dwpt,'Ty2'))  simple(diff(Dwpt,'dth1'))   simple(diff(Dwpt,'dph1'))  simple(diff(Dwpt,'dth2'))   simple(diff(Dwpt,'dph2'))  simple(diff(Dwpt,'th1'))   simple(diff(Dwpt,'ph1'))  simple(diff(Dwpt,'th2'))  simple(diff(Dwpt,'ph2'))  ];
disp('calculating F 8 of 12');
F8=[simple(diff(Dwyt,'Tp1'))   simple(diff(Dwyt,'Ty1'))  simple(diff(Dwyt,'Tp2'))   simple(diff(Dwyt,'Ty2'))  simple(diff(Dwyt,'dth1'))   simple(diff(Dwyt,'dph1'))  simple(diff(Dwyt,'dth2'))   simple(diff(Dwyt,'dph2'))  simple(diff(Dwyt,'th1'))   simple(diff(Dwyt,'ph1'))  simple(diff(Dwyt,'th2'))  simple(diff(Dwyt,'ph2'))  ];
disp('calculating F 9 of 12');
F9=[simple(diff(Dth1,'Tp1'))   simple(diff(Dth1,'Ty1'))  simple(diff(Dth1,'Tp2'))   simple(diff(Dth1,'Ty2'))  simple(diff(Dth1,'dth1'))   simple(diff(Dth1,'dph1'))  simple(diff(Dth1,'dth2'))   simple(diff(Dth1,'dph2'))  simple(diff(Dth1,'th1'))   simple(diff(Dth1,'ph1'))  simple(diff(Dth1,'th2'))  simple(diff(Dth1,'ph2'))  ];
disp('calculating F 10 of 12');
F10=[simple(diff(Dphi1,'Tp1'))  simple(diff(Dphi1,'Ty1')) simple(diff(Dphi1,'Tp2'))  simple(diff(Dphi1,'Ty2')) simple(diff(Dphi1,'dth1'))  simple(diff(Dphi1,'dph1')) simple(diff(Dphi1,'dth2'))  simple(diff(Dphi1,'dph2')) simple(diff(Dphi1,'th1'))  simple(diff(Dphi1,'ph1')) simple(diff(Dphi1,'th2')) simple(diff(Dphi1,'ph2')) ];
disp('calculating F 11 of 12');
F11=[simple(diff(Dth2,'Tp1'))   simple(diff(Dth2,'Ty1'))  simple(diff(Dth2,'Tp2'))   simple(diff(Dth2,'Ty2'))  simple(diff(Dth2,'dth1'))   simple(diff(Dth2,'dph1'))  simple(diff(Dth2,'dth2'))   simple(diff(Dth2,'dph2'))  simple(diff(Dth2,'th1'))   simple(diff(Dth2,'ph1'))  simple(diff(Dth2,'th2'))  simple(diff(Dth2,'ph2'))  ];
disp('calculating F 12 of 12');
F12=[simple(diff(Dphi2,'Tp1'))  simple(diff(Dphi2,'Ty1')) simple(diff(Dphi2,'Tp2'))  simple(diff(Dphi2,'Ty2')) simple(diff(Dphi2,'dth1'))  simple(diff(Dphi2,'dph1')) simple(diff(Dphi2,'dth2'))  simple(diff(Dphi2,'dph2')) simple(diff(Dphi2,'th1'))  simple(diff(Dphi2,'ph1')) simple(diff(Dphi2,'th2')) simple(diff(Dphi2,'ph2')) ];
F= [F1;F2;F3;F4;F5;F6;F7;F8;F9;F10;F11;F12]

