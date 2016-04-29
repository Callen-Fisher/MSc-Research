%code to plot the graphs (Cd vs speed plots)
figure(6)
plot(Fvel(:),CdBlackFur_0(:),'c+');
hold on
plot(Fvel(:),CdBlackFur_15(:),'k+');
plot(Fvel(:),CdBlackFur_30(:),'r+');
plot(Fvel(:),CdBlackFur_45(:),'m+');
plot(Fvel(:),CdBlackFur_60(:),'b+');
xlabel('speed (m/s)');
ylabel('Cd');
title('Cd vs speed for black fur at different angles of attack');
legend('black fur 0','black fur 15','black fur 30','black fur 45','black fur 60');

figure(7)
plot(Fvel(:),CdBlackFurSpray_0(:),'c+');
hold on
plot(Fvel(:),CdBlackFurSpray_15(:),'k+');
plot(Fvel(:),CdBlackFurSpray_30(:),'r+');
plot(Fvel(:),CdBlackFurSpray_45(:),'m+');
plot(Fvel(:),CdBlackFurSpray_60(:),'b+');
xlabel('speed (m/s)');
ylabel('Cd');
title('Cd vs speed for black fur with spray at different angles of attack');
legend('black fur spray 0','black fur spray 15','black fur spray 30','black fur spray 45','black fur spray 60');

figure(8)
plot(Fvel(:),CdCheetah_0(:),'c+');
hold on
plot(Fvel(:),CdCheetah_15(:),'k+');
plot(Fvel(:),CdCheetah_30(:),'r+');
plot(Fvel(:),CdCheetah_45(:),'m+');
plot(Fvel(:),CdCheetah_60(:),'b+');
xlabel('speed (m/s)');
ylabel('Cd');
title('Cd vs speed for cheetah fur at different angles of attack');
legend('cheetah fur 0','cheetah fur 15','cheetah fur 30','cheetah fur 45','cheetah fur 60');

figure(9)
plot(Fvel(:),CdCheetahSpray_0(:),'c+');
hold on
plot(Fvel(:),CdCheetahSpray_15(:),'k+');
plot(Fvel(:),CdCheetahSpray_30(:),'r+');
plot(Fvel(:),CdCheetahSpray_45(:),'m+');
plot(Fvel(:),CdCheetahSpray_60(:),'b+');
xlabel('speed (m/s)');
ylabel('Cd');
title('Cd vs speed for cheetah fur with spray at different angles of attack');
legend('cheetah fur Spray 0','cheetah fur Spray 15','cheetah fur Spray 30','cheetah fur Spray 45','cheetah fur Spray 60');








Fvel1=[10,15,20,25,30];
y1=[mean(CdBlackFur_0(1:9)),mean(CdBlackFur_0(10:18)),mean(CdBlackFur_0(19:27)),mean(CdBlackFur_0(28:36)),mean(CdBlackFur_0(37:45))];
y2=[mean(CdBlackFur_15(1:9)),mean(CdBlackFur_15(10:18)),mean(CdBlackFur_15(19:27)),mean(CdBlackFur_15(28:36)),mean(CdBlackFur_15(37:45))];
y3=[mean(CdBlackFur_30(1:9)),mean(CdBlackFur_30(10:18)),mean(CdBlackFur_30(19:27)),mean(CdBlackFur_30(28:36)),mean(CdBlackFur_30(37:45))];
y4=[mean(CdBlackFur_45(1:9)),mean(CdBlackFur_45(10:18)),mean(CdBlackFur_45(19:27)),mean(CdBlackFur_45(28:36)),mean(CdBlackFur_45(37:45))];
y5=[mean(CdBlackFur_60(1:9)),mean(CdBlackFur_60(10:18)),mean(CdBlackFur_60(19:27)),mean(CdBlackFur_60(28:36)),mean(CdBlackFur_60(37:45))];

u1=[max(CdBlackFur_0(1:9)),max(CdBlackFur_0(10:18)),max(CdBlackFur_0(19:27)),max(CdBlackFur_0(28:36)),max(CdBlackFur_0(37:45))]-y1;
u2=[max(CdBlackFur_15(1:9)),max(CdBlackFur_15(10:18)),max(CdBlackFur_15(19:27)),max(CdBlackFur_15(28:36)),max(CdBlackFur_15(37:45))]-y2;
u3=[max(CdBlackFur_30(1:9)),max(CdBlackFur_30(10:18)),max(CdBlackFur_30(19:27)),max(CdBlackFur_30(28:36)),max(CdBlackFur_30(37:45))]-y3;
u4=[max(CdBlackFur_45(1:9)),max(CdBlackFur_45(10:18)),max(CdBlackFur_45(19:27)),max(CdBlackFur_45(28:36)),max(CdBlackFur_45(37:45))]-y4;
u5=[max(CdBlackFur_60(1:9)),max(CdBlackFur_60(10:18)),max(CdBlackFur_60(19:27)),max(CdBlackFur_60(28:36)),max(CdBlackFur_60(37:45))]-y5;

l1=y1-[min(CdBlackFur_0(1:9)),min(CdBlackFur_0(10:18)),min(CdBlackFur_0(19:27)),min(CdBlackFur_0(28:36)),min(CdBlackFur_0(37:45))];
l2=y2-[min(CdBlackFur_15(1:9)),min(CdBlackFur_15(10:18)),min(CdBlackFur_15(19:27)),min(CdBlackFur_15(28:36)),min(CdBlackFur_15(37:45))];
l3=y3-[min(CdBlackFur_30(1:9)),min(CdBlackFur_30(10:18)),min(CdBlackFur_30(19:27)),min(CdBlackFur_30(28:36)),min(CdBlackFur_30(37:45))];
l4=y4-[min(CdBlackFur_45(1:9)),min(CdBlackFur_45(10:18)),min(CdBlackFur_45(19:27)),min(CdBlackFur_45(28:36)),min(CdBlackFur_45(37:45))];
l5=y5-[min(CdBlackFur_60(1:9)),min(CdBlackFur_60(10:18)),min(CdBlackFur_60(19:27)),min(CdBlackFur_60(28:36)),min(CdBlackFur_60(37:45))];

figure(20)
errorbar(Fvel1(:),y1(:),l1,u1,'c+');
hold on
errorbar(Fvel1(:),y2(:),l2,u2,'k+');
errorbar(Fvel1(:),y3(:),l3,u3,'r+');
errorbar(Fvel1(:),y4(:),l4,u4,'m+');
errorbar(Fvel1(:),y5(:),l5,u5,'b+');
plot(Fvel1(:),y1(:),'c');
plot(Fvel1(:),y2(:),'k');
plot(Fvel1(:),y3(:),'r');
plot(Fvel1(:),y4(:),'m');
plot(Fvel1(:),y5(:),'b');
xlabel('speed (m/s)');
ylabel('Cd');
title('Cd vs speed for black fur at different angles of attack');
legend('black fur 0','black fur 15','black fur 30','black fur 45','black fur 60');





Fvel1=[10,15,20,25,30];
y1=[mean(CdBlackFurSpray_0(1:9)),mean(CdBlackFurSpray_0(10:18)),mean(CdBlackFurSpray_0(19:27)),mean(CdBlackFurSpray_0(28:36)),mean(CdBlackFurSpray_0(37:45))];
y2=[mean(CdBlackFurSpray_15(1:9)),mean(CdBlackFurSpray_15(10:18)),mean(CdBlackFurSpray_15(19:27)),mean(CdBlackFurSpray_15(28:36)),mean(CdBlackFurSpray_15(37:45))];
y3=[mean(CdBlackFurSpray_30(1:9)),mean(CdBlackFurSpray_30(10:18)),mean(CdBlackFurSpray_30(19:27)),mean(CdBlackFurSpray_30(28:36)),mean(CdBlackFurSpray_30(37:45))];
y4=[mean(CdBlackFurSpray_45(1:9)),mean(CdBlackFurSpray_45(10:18)),mean(CdBlackFurSpray_45(19:27)),mean(CdBlackFurSpray_45(28:36)),mean(CdBlackFurSpray_45(37:45))];
y5=[mean(CdBlackFurSpray_60(1:9)),mean(CdBlackFurSpray_60(10:18)),mean(CdBlackFurSpray_60(19:27)),mean(CdBlackFurSpray_60(28:36)),mean(CdBlackFurSpray_60(37:45))];

u1=[max(CdBlackFurSpray_0(1:9)),max(CdBlackFurSpray_0(10:18)),max(CdBlackFurSpray_0(19:27)),max(CdBlackFurSpray_0(28:36)),max(CdBlackFurSpray_0(37:45))]-y1;
u2=[max(CdBlackFurSpray_15(1:9)),max(CdBlackFurSpray_15(10:18)),max(CdBlackFurSpray_15(19:27)),max(CdBlackFurSpray_15(28:36)),max(CdBlackFurSpray_15(37:45))]-y2;
u3=[max(CdBlackFurSpray_30(1:9)),max(CdBlackFurSpray_30(10:18)),max(CdBlackFurSpray_30(19:27)),max(CdBlackFurSpray_30(28:36)),max(CdBlackFurSpray_30(37:45))]-y3;
u4=[max(CdBlackFurSpray_45(1:9)),max(CdBlackFurSpray_45(10:18)),max(CdBlackFurSpray_45(19:27)),max(CdBlackFurSpray_45(28:36)),max(CdBlackFurSpray_45(37:45))]-y4;
u5=[max(CdBlackFurSpray_60(1:9)),max(CdBlackFurSpray_60(10:18)),max(CdBlackFurSpray_60(19:27)),max(CdBlackFurSpray_60(28:36)),max(CdBlackFurSpray_60(37:45))]-y5;

l1=y1-[min(CdBlackFurSpray_0(1:9)),min(CdBlackFurSpray_0(10:18)),min(CdBlackFurSpray_0(19:27)),min(CdBlackFurSpray_0(28:36)),min(CdBlackFurSpray_0(37:45))];
l2=y2-[min(CdBlackFurSpray_15(1:9)),min(CdBlackFurSpray_15(10:18)),min(CdBlackFurSpray_15(19:27)),min(CdBlackFurSpray_15(28:36)),min(CdBlackFurSpray_15(37:45))];
l3=y3-[min(CdBlackFurSpray_30(1:9)),min(CdBlackFurSpray_30(10:18)),min(CdBlackFurSpray_30(19:27)),min(CdBlackFurSpray_30(28:36)),min(CdBlackFurSpray_30(37:45))];
l4=y4-[min(CdBlackFurSpray_45(1:9)),min(CdBlackFurSpray_45(10:18)),min(CdBlackFurSpray_45(19:27)),min(CdBlackFurSpray_45(28:36)),min(CdBlackFurSpray_45(37:45))];
l5=y5-[min(CdBlackFurSpray_60(1:9)),min(CdBlackFurSpray_60(10:18)),min(CdBlackFurSpray_60(19:27)),min(CdBlackFurSpray_60(28:36)),min(CdBlackFurSpray_60(37:45))];




figure(21)
errorbar(Fvel1(:),y1(:),l1,u1,'c+');
hold on
errorbar(Fvel1(:),y2(:),l2,u2,'k+');
errorbar(Fvel1(:),y3(:),l3,u3,'r+');
errorbar(Fvel1(:),y4(:),l4,u4,'m+');
errorbar(Fvel1(:),y5(:),l5,u5,'b+');
plot(Fvel1(:),y1(:),'c');
plot(Fvel1(:),y2(:),'k');
plot(Fvel1(:),y3(:),'r');
plot(Fvel1(:),y4(:),'m');
plot(Fvel1(:),y5(:),'b');
xlabel('speed (m/s)');
ylabel('Cd');
title('Cd vs speed for black fur with spray at different angles of attack');
legend('black fur spray 0','black fur spray 15','black fur spray 30','black fur spray 45','black fur spray 60');


Fvel1=[10,15,20,25,30];
y1=[mean(CdCheetah_0(1:9)),mean(CdCheetah_0(10:18)),mean(CdCheetah_0(19:27)),mean(CdCheetah_0(28:36)),mean(CdCheetah_0(37:45))];
y2=[mean(CdCheetah_15(1:9)),mean(CdCheetah_15(10:18)),mean(CdCheetah_15(19:27)),mean(CdCheetah_15(28:36)),mean(CdCheetah_15(37:45))];
y3=[mean(CdCheetah_30(1:9)),mean(CdCheetah_30(10:18)),mean(CdCheetah_30(19:27)),mean(CdCheetah_30(28:36)),mean(CdCheetah_30(37:45))];
y4=[mean(CdCheetah_45(1:9)),mean(CdCheetah_45(10:18)),mean(CdCheetah_45(19:27)),mean(CdCheetah_45(28:36)),mean(CdCheetah_45(37:45))];
y5=[mean(CdCheetah_60(1:9)),mean(CdCheetah_60(10:18)),mean(CdCheetah_60(19:27)),mean(CdCheetah_60(28:36)),mean(CdCheetah_60(37:45))];

u1=[max(CdCheetah_0(1:9)),max(CdCheetah_0(10:18)),max(CdCheetah_0(19:27)),max(CdCheetah_0(28:36)),max(CdCheetah_0(37:45))]-y1;
u2=[max(CdCheetah_15(1:9)),max(CdCheetah_15(10:18)),max(CdCheetah_15(19:27)),max(CdCheetah_15(28:36)),max(CdCheetah_15(37:45))]-y2;
u3=[max(CdCheetah_30(1:9)),max(CdCheetah_30(10:18)),max(CdCheetah_30(19:27)),max(CdCheetah_30(28:36)),max(CdCheetah_30(37:45))]-y3;
u4=[max(CdCheetah_45(1:9)),max(CdCheetah_45(10:18)),max(CdCheetah_45(19:27)),max(CdCheetah_45(28:36)),max(CdCheetah_45(37:45))]-y4;
u5=[max(CdCheetah_60(1:9)),max(CdCheetah_60(10:18)),max(CdCheetah_60(19:27)),max(CdCheetah_60(28:36)),max(CdCheetah_60(37:45))]-y5;

l1=y1-[min(CdCheetah_0(1:9)),min(CdCheetah_0(10:18)),min(CdCheetah_0(19:27)),min(CdCheetah_0(28:36)),min(CdCheetah_0(37:45))];
l2=y2-[min(CdCheetah_15(1:9)),min(CdCheetah_15(10:18)),min(CdCheetah_15(19:27)),min(CdCheetah_15(28:36)),min(CdCheetah_15(37:45))];
l3=y3-[min(CdCheetah_30(1:9)),min(CdCheetah_30(10:18)),min(CdCheetah_30(19:27)),min(CdCheetah_30(28:36)),min(CdCheetah_30(37:45))];
l4=y4-[min(CdCheetah_45(1:9)),min(CdCheetah_45(10:18)),min(CdCheetah_45(19:27)),min(CdCheetah_45(28:36)),min(CdCheetah_45(37:45))];
l5=y5-[min(CdCheetah_60(1:9)),min(CdCheetah_60(10:18)),min(CdCheetah_60(19:27)),min(CdCheetah_60(28:36)),min(CdCheetah_60(37:45))];




figure(22)
errorbar(Fvel1(:),y1(:),l1,u1,'c+');
hold on
errorbar(Fvel1(:),y2(:),l2,u2,'k+');
errorbar(Fvel1(:),y3(:),l3,u3,'r+');
errorbar(Fvel1(:),y4(:),l4,u4,'m+');
errorbar(Fvel1(:),y5(:),l5,u5,'b+');
plot(Fvel1(:),y1(:),'c');
plot(Fvel1(:),y2(:),'k');
plot(Fvel1(:),y3(:),'r');
plot(Fvel1(:),y4(:),'m');
plot(Fvel1(:),y5(:),'b');
xlabel('speed (m/s)');
ylabel('Cd');
title('Cd vs speed for cheetah fur at different angles of attack');
legend('cheetah fur 0','cheetah fur 15','cheetah fur 30','cheetah fur 45','cheetah fur 60');




Fvel1=[10,15,20,25,30];
y1=[mean(CdCheetahSpray_0(1:9)),mean(CdCheetahSpray_0(10:18)),mean(CdCheetahSpray_0(19:27)),mean(CdCheetahSpray_0(28:36)),mean(CdCheetahSpray_0(37:45))];
y2=[mean(CdCheetahSpray_15(1:9)),mean(CdCheetahSpray_15(10:18)),mean(CdCheetahSpray_15(19:27)),mean(CdCheetahSpray_15(28:36)),mean(CdCheetahSpray_15(37:45))];
y3=[mean(CdCheetahSpray_30(1:9)),mean(CdCheetahSpray_30(10:18)),mean(CdCheetahSpray_30(19:27)),mean(CdCheetahSpray_30(28:36)),mean(CdCheetahSpray_30(37:45))];
y4=[mean(CdCheetahSpray_45(1:9)),mean(CdCheetahSpray_45(10:18)),mean(CdCheetahSpray_45(19:27)),mean(CdCheetahSpray_45(28:36)),mean(CdCheetahSpray_45(37:45))];
y5=[mean(CdCheetahSpray_60(1:9)),mean(CdCheetahSpray_60(10:18)),mean(CdCheetahSpray_60(19:27)),mean(CdCheetahSpray_60(28:36)),mean(CdCheetahSpray_60(37:45))];

u1=[max(CdCheetahSpray_0(1:9)),max(CdCheetahSpray_0(10:18)),max(CdCheetahSpray_0(19:27)),max(CdCheetahSpray_0(28:36)),max(CdCheetahSpray_0(37:45))]-y1;
u2=[max(CdCheetahSpray_15(1:9)),max(CdCheetahSpray_15(10:18)),max(CdCheetahSpray_15(19:27)),max(CdCheetahSpray_15(28:36)),max(CdCheetahSpray_15(37:45))]-y2;
u3=[max(CdCheetahSpray_30(1:9)),max(CdCheetahSpray_30(10:18)),max(CdCheetahSpray_30(19:27)),max(CdCheetahSpray_30(28:36)),max(CdCheetahSpray_30(37:45))]-y3;
u4=[max(CdCheetahSpray_45(1:9)),max(CdCheetahSpray_45(10:18)),max(CdCheetahSpray_45(19:27)),max(CdCheetahSpray_45(28:36)),max(CdCheetahSpray_45(37:45))]-y4;
u5=[max(CdCheetahSpray_60(1:9)),max(CdCheetahSpray_60(10:18)),max(CdCheetahSpray_60(19:27)),max(CdCheetahSpray_60(28:36)),max(CdCheetahSpray_60(37:45))]-y5;

l1=y1-[min(CdCheetahSpray_0(1:9)),min(CdCheetahSpray_0(10:18)),min(CdCheetahSpray_0(19:27)),min(CdCheetahSpray_0(28:36)),min(CdCheetahSpray_0(37:45))];
l2=y2-[min(CdCheetahSpray_15(1:9)),min(CdCheetahSpray_15(10:18)),min(CdCheetahSpray_15(19:27)),min(CdCheetahSpray_15(28:36)),min(CdCheetahSpray_15(37:45))];
l3=y3-[min(CdCheetahSpray_30(1:9)),min(CdCheetahSpray_30(10:18)),min(CdCheetahSpray_30(19:27)),min(CdCheetahSpray_30(28:36)),min(CdCheetahSpray_30(37:45))];
l4=y4-[min(CdCheetahSpray_45(1:9)),min(CdCheetahSpray_45(10:18)),min(CdCheetahSpray_45(19:27)),min(CdCheetahSpray_45(28:36)),min(CdCheetahSpray_45(37:45))];
l5=y5-[min(CdCheetahSpray_60(1:9)),min(CdCheetahSpray_60(10:18)),min(CdCheetahSpray_60(19:27)),min(CdCheetahSpray_60(28:36)),min(CdCheetahSpray_60(37:45))];






figure(23)
errorbar(Fvel1(:),y1(:),l1,u1,'c+');
hold on
errorbar(Fvel1(:),y2(:),l2,u2,'k+');
errorbar(Fvel1(:),y3(:),l3,u3,'r+');
errorbar(Fvel1(:),y4(:),l4,u4,'m+');
errorbar(Fvel1(:),y5(:),l5,u5,'b+');
plot(Fvel1(:),y1(:),'c');
plot(Fvel1(:),y2(:),'k');
plot(Fvel1(:),y3(:),'r');
plot(Fvel1(:),y4(:),'m');
plot(Fvel1(:),y5(:),'b');
xlabel('speed (m/s)');
ylabel('Cd');
title('Cd vs speed for cheetah fur with spray at different angles of attack');
legend('cheetah fur Spray 0','cheetah fur Spray 15','cheetah fur Spray 30','cheetah fur Spray 45','cheetah fur Spray 60');
