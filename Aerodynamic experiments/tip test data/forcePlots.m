%code to plot the force vs speed graphs
y1=[mean(FblackFur_0(1:9)),mean(FblackFur_0(10:18)),mean(FblackFur_0(19:27)),mean(FblackFur_0(28:36)),mean(FblackFur_0(37:45))];
y2=[mean(FblackFur_15(1:9)),mean(FblackFur_15(10:18)),mean(FblackFur_15(19:27)),mean(FblackFur_15(28:36)),mean(FblackFur_15(37:45))];
y3=[mean(FblackFur_30(1:9)),mean(FblackFur_30(10:18)),mean(FblackFur_30(19:27)),mean(FblackFur_30(28:36)),mean(FblackFur_30(37:45))];
y4=[mean(FblackFur_45(1:9)),mean(FblackFur_45(10:18)),mean(FblackFur_45(19:27)),mean(FblackFur_45(28:36)),mean(FblackFur_45(37:45))];
y5=[mean(FblackFur_60(1:9)),mean(FblackFur_60(10:18)),mean(FblackFur_60(19:27)),mean(FblackFur_60(28:36)),mean(FblackFur_60(37:45))];
Fvel1=[10,15,20,25,30];
figure(1)
plot(Fvel(:),FblackFur_0(:),'c+');
hold on
plot(Fvel(:),FblackFur_15(:),'k+');
plot(Fvel(:),FblackFur_30(:),'r+');
plot(Fvel(:),FblackFur_45(:),'m+');
plot(Fvel(:),FblackFur_60(:),'b+');
plot(Fvel1(:),y1(:),'c');
plot(Fvel1(:),y2(:),'k');
plot(Fvel1(:),y3(:),'r');
plot(Fvel1(:),y4(:),'m');
plot(Fvel1(:),y5(:),'b');
xlabel('speed (m/s)');
ylabel('F');
title('F vs speed for black fur at different angles of attack');
legend('black fur 0','black fur 15','black fur 30','black fur 45','black fur 60');


y1=[mean(FblackFurSpray_0(1:9)),mean(FblackFurSpray_0(10:18)),mean(FblackFurSpray_0(19:27)),mean(FblackFurSpray_0(28:36)),mean(FblackFurSpray_0(37:45))];
y2=[mean(FblackFurSpray_15(1:9)),mean(FblackFurSpray_15(10:18)),mean(FblackFurSpray_15(19:27)),mean(FblackFurSpray_15(28:36)),mean(FblackFurSpray_15(37:45))];
y3=[mean(FblackFurSpray_30(1:9)),mean(FblackFurSpray_30(10:18)),mean(FblackFurSpray_30(19:27)),mean(FblackFurSpray_30(28:36)),mean(FblackFurSpray_30(37:45))];
y4=[mean(FblackFurSpray_45(1:9)),mean(FblackFurSpray_45(10:18)),mean(FblackFurSpray_45(19:27)),mean(FblackFurSpray_45(28:36)),mean(FblackFurSpray_45(37:45))];
y5=[mean(FblackFurSpray_60(1:9)),mean(FblackFurSpray_60(10:18)),mean(FblackFurSpray_60(19:27)),mean(FblackFurSpray_60(28:36)),mean(FblackFurSpray_60(37:45))];


figure(2)
plot(Fvel(:),FblackFurSpray_0(:),'c+');
hold on
plot(Fvel(:),FblackFurSpray_15(:),'k+');
plot(Fvel(:),FblackFurSpray_30(:),'r+');
plot(Fvel(:),FblackFurSpray_45(:),'m+');
plot(Fvel(:),FblackFurSpray_60(:),'b+');
plot(Fvel1(:),y1(:),'c');
plot(Fvel1(:),y2(:),'k');
plot(Fvel1(:),y3(:),'r');
plot(Fvel1(:),y4(:),'m');
plot(Fvel1(:),y5(:),'b');
xlabel('speed (m/s)');
ylabel('F');
title('F vs speed for black fur with spray at different angles of attack');
legend('black fur spray 0','black fur spray 15','black fur spray 30','black fur spray 45','black fur spray 60');



y1=[mean(Fcheetah_0(1:9)),mean(Fcheetah_0(10:18)),mean(Fcheetah_0(19:27)),mean(Fcheetah_0(28:36)),mean(Fcheetah_0(37:45))];
y2=[mean(Fcheetah_15(1:9)),mean(Fcheetah_15(10:18)),mean(Fcheetah_15(19:27)),mean(Fcheetah_15(28:36)),mean(Fcheetah_15(37:45))];
y3=[mean(Fcheetah_30(1:9)),mean(Fcheetah_30(10:18)),mean(Fcheetah_30(19:27)),mean(Fcheetah_30(28:36)),mean(Fcheetah_30(37:45))];
y4=[mean(Fcheetah_45(1:9)),mean(Fcheetah_45(10:18)),mean(Fcheetah_45(19:27)),mean(Fcheetah_45(28:36)),mean(Fcheetah_45(37:45))];
y5=[mean(Fcheetah_60(1:9)),mean(Fcheetah_60(10:18)),mean(Fcheetah_60(19:27)),mean(Fcheetah_60(28:36)),mean(Fcheetah_60(37:45))];





figure(3)
plot(Fvel(:),Fcheetah_0(:),'c+');
hold on
plot(Fvel(:),Fcheetah_15(:),'k+');
plot(Fvel(:),Fcheetah_30(:),'r+');
plot(Fvel(:),Fcheetah_45(:),'m+');
plot(Fvel(:),Fcheetah_60(:),'b+');
plot(Fvel1(:),y1(:),'c');
plot(Fvel1(:),y2(:),'k');
plot(Fvel1(:),y3(:),'r');
plot(Fvel1(:),y4(:),'m');
plot(Fvel1(:),y5(:),'b');
xlabel('speed (m/s)');
ylabel('F');
title('F vs speed for cheetah fur at different angles of attack');
legend('cheetah fur 0','cheetah fur 15','cheetah fur 30','cheetah fur 45','cheetah fur 60');




y1=[mean(FcheetahSpray_0(1:9)),mean(FcheetahSpray_0(10:18)),mean(FcheetahSpray_0(19:27)),mean(FcheetahSpray_0(28:36)),mean(FcheetahSpray_0(37:45))];
y2=[mean(FcheetahSpray_15(1:9)),mean(FcheetahSpray_15(10:18)),mean(FcheetahSpray_15(19:27)),mean(FcheetahSpray_15(28:36)),mean(FcheetahSpray_15(37:45))];
y3=[mean(FcheetahSpray_30(1:9)),mean(FcheetahSpray_30(10:18)),mean(FcheetahSpray_30(19:27)),mean(FcheetahSpray_30(28:36)),mean(FcheetahSpray_30(37:45))];
y4=[mean(FcheetahSpray_45(1:9)),mean(FcheetahSpray_45(10:18)),mean(FcheetahSpray_45(19:27)),mean(FcheetahSpray_45(28:36)),mean(FcheetahSpray_45(37:45))];
y5=[mean(FcheetahSpray_60(1:9)),mean(FcheetahSpray_60(10:18)),mean(FcheetahSpray_60(19:27)),mean(FcheetahSpray_60(28:36)),mean(FcheetahSpray_60(37:45))];




figure(4)
plot(Fvel(:),FcheetahSpray_0(:),'c+');
hold on
plot(Fvel(:),FcheetahSpray_15(:),'k+');
plot(Fvel(:),FcheetahSpray_30(:),'r+');
plot(Fvel(:),FcheetahSpray_45(:),'m+');
plot(Fvel(:),FcheetahSpray_60(:),'b+');
plot(Fvel1(:),y1(:),'c');
plot(Fvel1(:),y2(:),'k');
plot(Fvel1(:),y3(:),'r');
plot(Fvel1(:),y4(:),'m');
plot(Fvel1(:),y5(:),'b');
xlabel('speed (m/s)');
ylabel('F');
title('F vs speed for cheetah fur with spray at different angles of attack');
legend('cheetah fur Spray 0','cheetah fur Spray 15','cheetah fur Spray 30','cheetah fur Spray 45','cheetah fur Spray 60');