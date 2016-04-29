%code used to generate the Cd vs reynolds plots 
%parameters
diameterFur=18.5*10^(-3);
diameterPipe=15*10^(-3);
diameterCheetah=16*10^(-3);
dynamicViscosity1=1.8124*10^(-5);
%reynolds numbers 
r1=[10,15,20,25,30]*rho*diameterFur/dynamicViscosity1;%fur
r2=[10,15,20,25,30]*rho*diameterPipe/dynamicViscosity1;%wood
r3=[10,15,20,25,30]*rho*diameterCheetah/dynamicViscosity1;%cheetah




figure(50)
y1=[mean(CdBlackFur_0(1:9)),mean(CdBlackFur_0(10:18)),mean(CdBlackFur_0(19:27)),mean(CdBlackFur_0(28:36)),mean(CdBlackFur_0(37:45))];
y2=[mean(CdBlackFurSpray_0(1:9)),mean(CdBlackFurSpray_0(10:18)),mean(CdBlackFurSpray_0(19:27)),mean(CdBlackFurSpray_0(28:36)),mean(CdBlackFurSpray_0(37:45))];%average Cd for fur 
y4=[mean(CdCheetah_0(1:9)),mean(CdCheetah_0(10:18)),mean(CdCheetah_0(19:27)),mean(CdCheetah_0(28:36)),mean(CdCheetah_0(37:45))];
y5=[mean(CdCheetahSpray_0(1:9)),mean(CdCheetahSpray_0(10:18)),mean(CdCheetahSpray_0(19:27)),mean(CdCheetahSpray_0(28:36)),mean(CdCheetahSpray_0(37:45))];
loglog(r1,y1,'c');
hold on
loglog(r1,y2,'k');
%loglog(r2,y3,'r');
loglog(r3,y4,'m');
loglog(r3,y5,'b');
xlabel('Reynolds number');
ylabel('Cd');
title('Cd vs reynolds number for 0 degrees');
legend('fur','spray','cheetah','cheetah spray');

figure(51)
y1=[mean(CdBlackFur_15(1:9)),mean(CdBlackFur_15(10:18)),mean(CdBlackFur_15(19:27)),mean(CdBlackFur_15(28:36)),mean(CdBlackFur_15(37:45))];
y2=[mean(CdBlackFurSpray_15(1:9)),mean(CdBlackFurSpray_15(10:18)),mean(CdBlackFurSpray_15(19:27)),mean(CdBlackFurSpray_15(28:36)),mean(CdBlackFurSpray_15(37:45))];
y4=[mean(CdCheetah_15(1:9)),mean(CdCheetah_15(10:18)),mean(CdCheetah_15(19:27)),mean(CdCheetah_15(28:36)),mean(CdCheetah_15(37:45))];
y5=[mean(CdCheetahSpray_15(1:9)),mean(CdCheetahSpray_15(10:18)),mean(CdCheetahSpray_15(19:27)),mean(CdCheetahSpray_15(28:36)),mean(CdCheetahSpray_15(37:45))];
loglog(r1,y1,'c');
hold on
loglog(r1,y2,'k');
loglog(r3,y4,'m');
loglog(r3,y5,'b');
xlabel('Reynolds number');
ylabel('Cd');
title('Cd vs reynolds number for 15 degrees');
legend('fur','spray','cheetah','cheetah spray');

figure(52)
y1=[mean(CdBlackFur_30(1:9)),mean(CdBlackFur_30(10:18)),mean(CdBlackFur_30(19:27)),mean(CdBlackFur_30(28:36)),mean(CdBlackFur_30(37:45))];
y2=[mean(CdBlackFurSpray_30(1:9)),mean(CdBlackFurSpray_30(10:18)),mean(CdBlackFurSpray_30(19:27)),mean(CdBlackFurSpray_30(28:36)),mean(CdBlackFurSpray_30(37:45))];
y4=[mean(CdCheetah_30(1:9)),mean(CdCheetah_30(10:18)),mean(CdCheetah_30(19:27)),mean(CdCheetah_30(28:36)),mean(CdCheetah_30(37:45))];
y5=[mean(CdCheetahSpray_30(1:9)),mean(CdCheetahSpray_30(10:18)),mean(CdCheetahSpray_30(19:27)),mean(CdCheetahSpray_30(28:36)),mean(CdCheetahSpray_30(37:45))];
loglog(r1,y1,'c');
hold on
loglog(r1,y2,'k');
loglog(r3,y4,'m');
loglog(r3,y5,'b');
xlabel('Reynolds number');
ylabel('Cd');
title('Cd vs reynolds number for 30 degrees');
legend('fur','spray','cheetah','cheetah spray');

figure(53)
y1=[mean(CdBlackFur_45(1:9)),mean(CdBlackFur_45(10:18)),mean(CdBlackFur_45(19:27)),mean(CdBlackFur_45(28:36)),mean(CdBlackFur_45(37:45))];
y2=[mean(CdBlackFurSpray_45(1:9)),mean(CdBlackFurSpray_45(10:18)),mean(CdBlackFurSpray_45(19:27)),mean(CdBlackFurSpray_45(28:36)),mean(CdBlackFurSpray_45(37:45))];
y4=[mean(CdCheetah_45(1:9)),mean(CdCheetah_45(10:18)),mean(CdCheetah_45(19:27)),mean(CdCheetah_45(28:36)),mean(CdCheetah_45(37:45))];
y5=[mean(CdCheetahSpray_45(1:9)),mean(CdCheetahSpray_45(10:18)),mean(CdCheetahSpray_45(19:27)),mean(CdCheetahSpray_45(28:36)),mean(CdCheetahSpray_45(37:45))];
loglog(r1,y1,'c');
hold on
loglog(r1,y2,'k');
loglog(r3,y4,'m');
loglog(r3,y5,'b');
xlabel('Reynolds number');
ylabel('Cd');
title('Cd vs reynolds number for 45 degrees');
legend('fur','spray','cheetah','cheetah spray');

figure(54)
y1=[mean(CdBlackFur_60(1:9)),mean(CdBlackFur_60(10:18)),mean(CdBlackFur_60(19:27)),mean(CdBlackFur_60(28:36)),mean(CdBlackFur_60(37:45))];
y2=[mean(CdBlackFurSpray_60(1:9)),mean(CdBlackFurSpray_60(10:18)),mean(CdBlackFurSpray_60(19:27)),mean(CdBlackFurSpray_60(28:36)),mean(CdBlackFurSpray_60(37:45))];
y4=[mean(CdCheetah_60(1:9)),mean(CdCheetah_60(10:18)),mean(CdCheetah_60(19:27)),mean(CdCheetah_60(28:36)),mean(CdCheetah_60(37:45))];
y5=[mean(CdCheetahSpray_60(1:9)),mean(CdCheetahSpray_60(10:18)),mean(CdCheetahSpray_60(19:27)),mean(CdCheetahSpray_60(28:36)),mean(CdCheetahSpray_60(37:45))];
loglog(r1,y1,'c');
hold on
loglog(r1,y2,'k');
loglog(r3,y4,'m');
loglog(r3,y5,'b');
xlabel('Reynolds number');
ylabel('Cd');
title('Cd vs reynolds number for 60 degrees');
legend('fur','spray','cheetah','cheetah spray');