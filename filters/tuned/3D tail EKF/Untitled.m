l1=2;
th1=1.1157;
ph1=-0.1235;

Rpitch1=[[cos(th1) 0 sin(th1)];[0 1 0];[-sin(th1) 0 cos(th1)]];
Ryaw1=[[cos(ph1) -sin(ph1) 0];[sin(ph1) cos(ph1) 0];[0 0 1]];
R_0_1=(Ryaw1*Rpitch1);
P1=R_0_1*[l1/2; 0; 0]