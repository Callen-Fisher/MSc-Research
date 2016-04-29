%code to plot the figures (Cd vs angle)
figure(11)
plot(angle(:),CdBlackFur10(:),'c+');
hold on
plot(angle(:),CdBlackFur15(:),'k+');
plot(angle(:),CdBlackFur20(:),'r+');
plot(angle(:),CdBlackFur25(:),'m+');
plot(angle(:),CdBlackFur30(:),'b+');
xlabel('angle');
ylabel('Cd');
title('Cd vs angle for black fur');
legend('10m/s','15m/s','20m/s','25m/s','30m/s');

figure(12)
plot(angle(:),CdBlackFurSpray10(:),'c+');
hold on
plot(angle(:),CdBlackFurSpray15(:),'k+');
plot(angle(:),CdBlackFurSpray20(:),'r+');
plot(angle(:),CdBlackFurSpray25(:),'m+');
plot(angle(:),CdBlackFurSpray30(:),'b+');
xlabel('angle');
ylabel('Cd');
title('Cd vs angle for black fur with spray');
legend('10m/s','15m/s','20m/s','25m/s','30m/s');

figure(13)
plot(angle(:),CdCheetah10(:),'c+');
hold on
plot(angle(:),CdCheetah15(:),'k+');
plot(angle(:),CdCheetah20(:),'r+');
plot(angle(:),CdCheetah25(:),'m+');
plot(angle(:),CdCheetah30(:),'b+');
xlabel('angle');
ylabel('Cd');
title('Cd vs angle for cheetah fur');
legend('10m/s','15m/s','20m/s','25m/s','30m/s');

figure(14)
plot(angle(:),CdCheetahSpray10(:),'c+');
hold on
plot(angle(:),CdCheetahSpray15(:),'k+');
plot(angle(:),CdCheetahSpray20(:),'r+');
plot(angle(:),CdCheetahSpray25(:),'m+');
plot(angle(:),CdCheetahSpray30(:),'b+');
xlabel('angle');
ylabel('Cd');
title('Cd vs angle for cheetah fur with spray');
legend('10m/s','15m/s','20m/s','25m/s','30m/s');










angle1=[0,15,30,45,60];
y1=[mean(CdBlackFur10(1:9)),mean(CdBlackFur10(10:18)),mean(CdBlackFur10(19:27)),mean(CdBlackFur10(28:36)),mean(CdBlackFur10(37:45))];
y2=[mean(CdBlackFur15(1:9)),mean(CdBlackFur15(10:18)),mean(CdBlackFur15(19:27)),mean(CdBlackFur15(28:36)),mean(CdBlackFur15(37:45))];
y3=[mean(CdBlackFur20(1:9)),mean(CdBlackFur20(10:18)),mean(CdBlackFur20(19:27)),mean(CdBlackFur20(28:36)),mean(CdBlackFur20(37:45))];
y4=[mean(CdBlackFur25(1:9)),mean(CdBlackFur25(10:18)),mean(CdBlackFur25(19:27)),mean(CdBlackFur25(28:36)),mean(CdBlackFur25(37:45))];
y5=[mean(CdBlackFur30(1:9)),mean(CdBlackFur30(10:18)),mean(CdBlackFur30(19:27)),mean(CdBlackFur30(28:36)),mean(CdBlackFur30(37:45))];

u1=[max(CdBlackFur10(1:9)),max(CdBlackFur10(10:18)),max(CdBlackFur10(19:27)),max(CdBlackFur10(28:36)),max(CdBlackFur10(37:45))]-y1;
u2=[max(CdBlackFur15(1:9)),max(CdBlackFur15(10:18)),max(CdBlackFur15(19:27)),max(CdBlackFur15(28:36)),max(CdBlackFur15(37:45))]-y2;
u3=[max(CdBlackFur20(1:9)),max(CdBlackFur20(10:18)),max(CdBlackFur20(19:27)),max(CdBlackFur20(28:36)),max(CdBlackFur20(37:45))]-y3;
u4=[max(CdBlackFur25(1:9)),max(CdBlackFur25(10:18)),max(CdBlackFur25(19:27)),max(CdBlackFur25(28:36)),max(CdBlackFur25(37:45))]-y4;
u5=[max(CdBlackFur30(1:9)),max(CdBlackFur30(10:18)),max(CdBlackFur30(19:27)),max(CdBlackFur30(28:36)),max(CdBlackFur30(37:45))]-y5;

l1=y1-[min(CdBlackFur10(1:9)),min(CdBlackFur10(10:18)),min(CdBlackFur10(19:27)),min(CdBlackFur10(28:36)),min(CdBlackFur10(37:45))];
l2=y2-[min(CdBlackFur15(1:9)),min(CdBlackFur15(10:18)),min(CdBlackFur15(19:27)),min(CdBlackFur15(28:36)),min(CdBlackFur15(37:45))];
l3=y3-[min(CdBlackFur20(1:9)),min(CdBlackFur20(10:18)),min(CdBlackFur20(19:27)),min(CdBlackFur20(28:36)),min(CdBlackFur20(37:45))];
l4=y4-[min(CdBlackFur25(1:9)),min(CdBlackFur25(10:18)),min(CdBlackFur25(19:27)),min(CdBlackFur25(28:36)),min(CdBlackFur25(37:45))];
l5=y5-[min(CdBlackFur30(1:9)),min(CdBlackFur30(10:18)),min(CdBlackFur30(19:27)),min(CdBlackFur30(28:36)),min(CdBlackFur30(37:45))];



figure(30)
errorbar(angle1(:),y1(:),l1,u1,'c+');
hold on
errorbar(angle1(:),y2(:),l2,u2,'k+');
errorbar(angle1(:),y3(:),l3,u3,'r+');
errorbar(angle1(:),y4(:),l4,u4,'m+');
errorbar(angle1(:),y5(:),l5,u5,'b+');
plot(angle1(:),y1(:),'c');
plot(angle1(:),y2(:),'k');
plot(angle1(:),y3(:),'r');
plot(angle1(:),y4(:),'m');
plot(angle1(:),y5(:),'b');
xlabel('angle');
ylabel('Cd');
title('Cd vs angle for black fur');
legend('10m/s','15m/s','20m/s','25m/s','30m/s');





angle1=[0,15,30,45,60];
y1=[mean(CdBlackFurSpray10(1:9)),mean(CdBlackFurSpray10(10:18)),mean(CdBlackFurSpray10(19:27)),mean(CdBlackFurSpray10(28:36)),mean(CdBlackFurSpray10(37:45))];
y2=[mean(CdBlackFurSpray15(1:9)),mean(CdBlackFurSpray15(10:18)),mean(CdBlackFurSpray15(19:27)),mean(CdBlackFurSpray15(28:36)),mean(CdBlackFurSpray15(37:45))];
y3=[mean(CdBlackFurSpray20(1:9)),mean(CdBlackFurSpray20(10:18)),mean(CdBlackFurSpray20(19:27)),mean(CdBlackFurSpray20(28:36)),mean(CdBlackFurSpray20(37:45))];
y4=[mean(CdBlackFurSpray25(1:9)),mean(CdBlackFurSpray25(10:18)),mean(CdBlackFurSpray25(19:27)),mean(CdBlackFurSpray25(28:36)),mean(CdBlackFurSpray25(37:45))];
y5=[mean(CdBlackFurSpray30(1:9)),mean(CdBlackFurSpray30(10:18)),mean(CdBlackFurSpray30(19:27)),mean(CdBlackFurSpray30(28:36)),mean(CdBlackFurSpray30(37:45))];

u1=[max(CdBlackFurSpray10(1:9)),max(CdBlackFurSpray10(10:18)),max(CdBlackFurSpray10(19:27)),max(CdBlackFurSpray10(28:36)),max(CdBlackFurSpray10(37:45))]-y1;
u2=[max(CdBlackFurSpray15(1:9)),max(CdBlackFurSpray15(10:18)),max(CdBlackFurSpray15(19:27)),max(CdBlackFurSpray15(28:36)),max(CdBlackFurSpray15(37:45))]-y2;
u3=[max(CdBlackFurSpray20(1:9)),max(CdBlackFurSpray20(10:18)),max(CdBlackFurSpray20(19:27)),max(CdBlackFurSpray20(28:36)),max(CdBlackFurSpray20(37:45))]-y3;
u4=[max(CdBlackFurSpray25(1:9)),max(CdBlackFurSpray25(10:18)),max(CdBlackFurSpray25(19:27)),max(CdBlackFurSpray25(28:36)),max(CdBlackFurSpray25(37:45))]-y4;
u5=[max(CdBlackFurSpray30(1:9)),max(CdBlackFurSpray30(10:18)),max(CdBlackFurSpray30(19:27)),max(CdBlackFurSpray30(28:36)),max(CdBlackFurSpray30(37:45))]-y5;

l1=y1-[min(CdBlackFurSpray10(1:9)),min(CdBlackFurSpray10(10:18)),min(CdBlackFurSpray10(19:27)),min(CdBlackFurSpray10(28:36)),min(CdBlackFurSpray10(37:45))];
l2=y2-[min(CdBlackFurSpray15(1:9)),min(CdBlackFurSpray15(10:18)),min(CdBlackFurSpray15(19:27)),min(CdBlackFurSpray15(28:36)),min(CdBlackFurSpray15(37:45))];
l3=y3-[min(CdBlackFurSpray20(1:9)),min(CdBlackFurSpray20(10:18)),min(CdBlackFurSpray20(19:27)),min(CdBlackFurSpray20(28:36)),min(CdBlackFurSpray20(37:45))];
l4=y4-[min(CdBlackFurSpray25(1:9)),min(CdBlackFurSpray25(10:18)),min(CdBlackFurSpray25(19:27)),min(CdBlackFurSpray25(28:36)),min(CdBlackFurSpray25(37:45))];
l5=y5-[min(CdBlackFurSpray30(1:9)),min(CdBlackFurSpray30(10:18)),min(CdBlackFurSpray30(19:27)),min(CdBlackFurSpray30(28:36)),min(CdBlackFurSpray30(37:45))];






figure(31)
errorbar(angle1(:),y1(:),l1,u1,'c+');
hold on
errorbar(angle1(:),y2(:),l2,u2,'k+');
errorbar(angle1(:),y3(:),l3,u3,'r+');
errorbar(angle1(:),y4(:),l4,u4,'m+');
errorbar(angle1(:),y5(:),l5,u5,'b+');
plot(angle1(:),y1(:),'c');
plot(angle1(:),y2(:),'k');
plot(angle1(:),y3(:),'r');
plot(angle1(:),y4(:),'m');
plot(angle1(:),y5(:),'b');
xlabel('angle');
ylabel('Cd');
title('Cd vs angle for black fur with spray');
legend('10m/s','15m/s','20m/s','25m/s','30m/s');



angle1=[0,15,30,45,60];
y1=[mean(CdCheetah10(1:9)),mean(CdCheetah10(10:18)),mean(CdCheetah10(19:27)),mean(CdCheetah10(28:36)),mean(CdCheetah10(37:45))]
y2=[mean(CdCheetah15(1:9)),mean(CdCheetah15(10:18)),mean(CdCheetah15(19:27)),mean(CdCheetah15(28:36)),mean(CdCheetah15(37:45))]
y3=[mean(CdCheetah20(1:9)),mean(CdCheetah20(10:18)),mean(CdCheetah20(19:27)),mean(CdCheetah20(28:36)),mean(CdCheetah20(37:45))]
y4=[mean(CdCheetah25(1:9)),mean(CdCheetah25(10:18)),mean(CdCheetah25(19:27)),mean(CdCheetah25(28:36)),mean(CdCheetah25(37:45))]
y5=[mean(CdCheetah30(1:9)),mean(CdCheetah30(10:18)),mean(CdCheetah30(19:27)),mean(CdCheetah30(28:36)),mean(CdCheetah30(37:45))]

u1=[max(CdCheetah10(1:9)),max(CdCheetah10(10:18)),max(CdCheetah10(19:27)),max(CdCheetah10(28:36)),max(CdCheetah10(37:45))]-y1;
u2=[max(CdCheetah15(1:9)),max(CdCheetah15(10:18)),max(CdCheetah15(19:27)),max(CdCheetah15(28:36)),max(CdCheetah15(37:45))]-y2;
u3=[max(CdCheetah20(1:9)),max(CdCheetah20(10:18)),max(CdCheetah20(19:27)),max(CdCheetah20(28:36)),max(CdCheetah20(37:45))]-y3;
u4=[max(CdCheetah25(1:9)),max(CdCheetah25(10:18)),max(CdCheetah25(19:27)),max(CdCheetah25(28:36)),max(CdCheetah25(37:45))]-y4;
u5=[max(CdCheetah30(1:9)),max(CdCheetah30(10:18)),max(CdCheetah30(19:27)),max(CdCheetah30(28:36)),max(CdCheetah30(37:45))]-y5;

l1=y1-[min(CdCheetah10(1:9)),min(CdCheetah10(10:18)),min(CdCheetah10(19:27)),min(CdCheetah10(28:36)),min(CdCheetah10(37:45))];
l2=y2-[min(CdCheetah15(1:9)),min(CdCheetah15(10:18)),min(CdCheetah15(19:27)),min(CdCheetah15(28:36)),min(CdCheetah15(37:45))];
l3=y3-[min(CdCheetah20(1:9)),min(CdCheetah20(10:18)),min(CdCheetah20(19:27)),min(CdCheetah20(28:36)),min(CdCheetah20(37:45))];
l4=y4-[min(CdCheetah25(1:9)),min(CdCheetah25(10:18)),min(CdCheetah25(19:27)),min(CdCheetah25(28:36)),min(CdCheetah25(37:45))];
l5=y5-[min(CdCheetah30(1:9)),min(CdCheetah30(10:18)),min(CdCheetah30(19:27)),min(CdCheetah30(28:36)),min(CdCheetah30(37:45))];

for i=1:length(y1)
    y1(i)=y1(i)/cos(angle1(i)*pi/180);
    y2(i)=y2(i)/cos(angle1(i)*pi/180);
    y3(i)=y3(i)/cos(angle1(i)*pi/180);
    y4(i)=y4(i)/cos(angle1(i)*pi/180);
    y5(i)=y5(i)/cos(angle1(i)*pi/180);
end

figure1=figure(32)
% errorbar(angle1(:),y1(:),l1,u1,'c+');
hold on
% errorbar(angle1(:),y2(:),l2,u2,'k+');
% errorbar(angle1(:),y3(:),l3,u3,'r+');
% errorbar(angle1(:),y4(:),l4,u4,'m+');
% errorbar(angle1(:),y5(:),l5,u5,'b+');
plot(angle1(:),y1(:),'Color',[0.0392156876623631 0.141176477074623 0.415686279535294],'LineWidth',2);%y1/(4.054/1.064)
plot(angle1(:),y2(:),'Color',[0.0431372560560703 0.517647087574005 0.780392169952393],'LineWidth',2);
plot(angle1(:),y3(:),'Color',[0.847058832645416 0.160784319043159 0],'LineWidth',2);
plot(angle1(:),y4(:),'Color',[0 0.749019622802734 0.749019622802734],'LineWidth',2);
plot(angle1(:),y5(:),'Color',[0.313725501298904 0.313725501298904 0.313725501298904],'LineWidth',2);
xlabel('Angle (\circ)','FontSize',16);
ylabel('Drag Coefficient (Cd)','FontSize',14);
title('Drag Coefficient vs angle of attack for cheetah fur','FontSize',14);
legend('v = 10m/s','v = 15m/s','v = 20m/s','v = 25m/s','v = 30m/s');
box on
grid on
figure(200)
yy=[y1(:), y2(:),y3(:),y4(:),y5(:)];

plot([10,15,20,25,30],yy')

angle1=[0,15,30,45,60];
y1=[mean(CdCheetahSpray10(1:9)),mean(CdCheetahSpray10(10:18)),mean(CdCheetahSpray10(19:27)),mean(CdCheetahSpray10(28:36)),mean(CdCheetahSpray10(37:45))];
y2=[mean(CdCheetahSpray15(1:9)),mean(CdCheetahSpray15(10:18)),mean(CdCheetahSpray15(19:27)),mean(CdCheetahSpray15(28:36)),mean(CdCheetahSpray15(37:45))];
y3=[mean(CdCheetahSpray20(1:9)),mean(CdCheetahSpray20(10:18)),mean(CdCheetahSpray20(19:27)),mean(CdCheetahSpray20(28:36)),mean(CdCheetahSpray20(37:45))];
y4=[mean(CdCheetahSpray25(1:9)),mean(CdCheetahSpray25(10:18)),mean(CdCheetahSpray25(19:27)),mean(CdCheetahSpray25(28:36)),mean(CdCheetahSpray25(37:45))];
y5=[mean(CdCheetahSpray30(1:9)),mean(CdCheetahSpray30(10:18)),mean(CdCheetahSpray30(19:27)),mean(CdCheetahSpray30(28:36)),mean(CdCheetahSpray30(37:45))];

u1=[max(CdCheetahSpray10(1:9)),max(CdCheetahSpray10(10:18)),max(CdCheetahSpray10(19:27)),max(CdCheetahSpray10(28:36)),max(CdCheetahSpray10(37:45))]-y1;
u2=[max(CdCheetahSpray15(1:9)),max(CdCheetahSpray15(10:18)),max(CdCheetahSpray15(19:27)),max(CdCheetahSpray15(28:36)),max(CdCheetahSpray15(37:45))]-y2;
u3=[max(CdCheetahSpray20(1:9)),max(CdCheetahSpray20(10:18)),max(CdCheetahSpray20(19:27)),max(CdCheetahSpray20(28:36)),max(CdCheetahSpray20(37:45))]-y3;
u4=[max(CdCheetahSpray25(1:9)),max(CdCheetahSpray25(10:18)),max(CdCheetahSpray25(19:27)),max(CdCheetahSpray25(28:36)),max(CdCheetahSpray25(37:45))]-y4;
u5=[max(CdCheetahSpray30(1:9)),max(CdCheetahSpray30(10:18)),max(CdCheetahSpray30(19:27)),max(CdCheetahSpray30(28:36)),max(CdCheetahSpray30(37:45))]-y5;

l1=y1-[min(CdCheetahSpray10(1:9)),min(CdCheetahSpray10(10:18)),min(CdCheetahSpray10(19:27)),min(CdCheetahSpray10(28:36)),min(CdCheetahSpray10(37:45))];
l2=y2-[min(CdCheetahSpray15(1:9)),min(CdCheetahSpray15(10:18)),min(CdCheetahSpray15(19:27)),min(CdCheetahSpray15(28:36)),min(CdCheetahSpray15(37:45))];
l3=y3-[min(CdCheetahSpray20(1:9)),min(CdCheetahSpray20(10:18)),min(CdCheetahSpray20(19:27)),min(CdCheetahSpray20(28:36)),min(CdCheetahSpray20(37:45))];
l4=y4-[min(CdCheetahSpray25(1:9)),min(CdCheetahSpray25(10:18)),min(CdCheetahSpray25(19:27)),min(CdCheetahSpray25(28:36)),min(CdCheetahSpray25(37:45))];
l5=y5-[min(CdCheetahSpray30(1:9)),min(CdCheetahSpray30(10:18)),min(CdCheetahSpray30(19:27)),min(CdCheetahSpray30(28:36)),min(CdCheetahSpray30(37:45))];



figure(33)
errorbar(angle1(:),y1(:),l1,u1,'c+');
hold on
errorbar(angle1(:),y2(:),l2,u2,'k+');
errorbar(angle1(:),y3(:),l3,u3,'r+');
errorbar(angle1(:),y4(:),l4,u4,'m+');
errorbar(angle1(:),y5(:),l5,u5,'b+');
plot(angle1(:),y1(:),'c');
plot(angle1(:),y2(:),'k');
plot(angle1(:),y3(:),'r');
plot(angle1(:),y4(:),'m');
plot(angle1(:),y5(:),'b');
xlabel('angle');
ylabel('Cd');
title('Cd vs angle for cheetah fur with spray');
legend('10m/s','15m/s','20m/s','25m/s','30m/s');


























angle1=[0,15,30,45,60];
y1=[mean(CdWood10(1:9)),mean(CdWood10(10:18)),mean(CdWood10(19:27)),mean(CdWood10(28:36)),mean(CdWood10(37:45))];
y2=[mean(CdWood15(1:9)),mean(CdWood15(10:18)),mean(CdWood15(19:27)),mean(CdWood15(28:36)),mean(CdWood15(37:45))];
y3=[mean(CdWood20(1:9)),mean(CdWood20(10:18)),mean(CdWood20(19:27)),mean(CdWood20(28:36)),mean(CdWood20(37:45))];
y4=[mean(CdWood25(1:9)),mean(CdWood25(10:18)),mean(CdWood25(19:27)),mean(CdWood25(28:36)),mean(CdWood25(37:45))];
y5=[mean(CdWood30(1:9)),mean(CdWood30(10:18)),mean(CdWood30(19:27)),mean(CdWood30(28:36)),mean(CdWood30(37:45))];

u1=[max(CdWood10(1:9)),max(CdWood10(10:18)),max(CdWood10(19:27)),max(CdWood10(28:36)),max(CdWood10(37:45))]-y1;
u2=[max(CdWood15(1:9)),max(CdWood15(10:18)),max(CdWood15(19:27)),max(CdWood15(28:36)),max(CdWood15(37:45))]-y2;
u3=[max(CdWood20(1:9)),max(CdWood20(10:18)),max(CdWood20(19:27)),max(CdWood20(28:36)),max(CdWood20(37:45))]-y3;
u4=[max(CdWood25(1:9)),max(CdWood25(10:18)),max(CdWood25(19:27)),max(CdWood25(28:36)),max(CdWood25(37:45))]-y4;
u5=[max(CdWood30(1:9)),max(CdWood30(10:18)),max(CdWood30(19:27)),max(CdWood30(28:36)),max(CdWood30(37:45))]-y5;

l1=y1-[min(CdWood10(1:9)),min(CdWood10(10:18)),min(CdWood10(19:27)),min(CdWood10(28:36)),min(CdWood10(37:45))];
l2=y2-[min(CdWood15(1:9)),min(CdWood15(10:18)),min(CdWood15(19:27)),min(CdWood15(28:36)),min(CdWood15(37:45))];
l3=y3-[min(CdWood20(1:9)),min(CdWood20(10:18)),min(CdWood20(19:27)),min(CdWood20(28:36)),min(CdWood20(37:45))];
l4=y4-[min(CdWood25(1:9)),min(CdWood25(10:18)),min(CdWood25(19:27)),min(CdWood25(28:36)),min(CdWood25(37:45))];
l5=y5-[min(CdWood30(1:9)),min(CdWood30(10:18)),min(CdWood30(19:27)),min(CdWood30(28:36)),min(CdWood30(37:45))];










figure(100)
errorbar(angle1(:),y1(:),l1,u1,'c+');
hold on
errorbar(angle1(:),y2(:),l2,u2,'k+');
errorbar(angle1(:),y3(:),l3,u3,'r+');
errorbar(angle1(:),y4(:),l4,u4,'m+');
errorbar(angle1(:),y5(:),l5,u5,'b+');
plot(angle1(:),y1(:),'c');
plot(angle1(:),y2(:),'k');
plot(angle1(:),y3(:),'r');
plot(angle1(:),y4(:),'m');
plot(angle1(:),y5(:),'b');
xlabel('angle');
ylabel('Cd');
title('Cd vs angle for wood');
legend('10m/s','15m/s','20m/s','25m/s','30m/s');