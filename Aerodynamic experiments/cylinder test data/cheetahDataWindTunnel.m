%cheetah fur- 8may cylinder test-19 degrees with pressure 1010.16
%other fur- 7may cylinder test with 17.2 degrees
%clc
%clear
%import the experimental data 
blackFur10=[.305,.303,.302,.305,.302,.305,.302,.304,.303];
blackFur15=[.615,.618,.619,.622,.616,.613,.615,.609,.611];
blackFur20=[1.042,1.037,1.055,1.046,1.055,1.047,1.046,1.049,1.064];
blackFur25=[1.401,1.415,1.414,1.411,1.393,1.405,1.402,1.402,1.405];
blackFur30=[1.840,1.833,1.841,1.847,1.842,1.831,1.864,1.851,1.831];

blackSpray10=[.317,.316,.319,.318,.316,.318,.316,.319,.317];
blackSpray15=[.646,.643,.641,.641,.652,.640,.648,.650,.643];
blackSpray20=[.939,.937,.935,.931,.939,.936,.940,.935,.940];
blackSpray25=[1.228,1.226,1.224,1.225,1.229,1.224,1.230,1.228,1.231];
blackSpray30=[1.512,1.516,1.513,1.514,1.521,1.495,1.510,1.494,1.504];

pipe10=[.137,.138,.135,.136,.139,.136,.139,.140,.138];
pipe15=[.318,.317,.323,.323,.322,.317,.319,.314,.326];
pipe20=[.545,.549,.551,.562,.557,.570,.561,.552,.553];
pipe25=[.721,.725,.750,.736,.739,.733,.748,.751,.750];
pipe30=[1.014,1.015,1.020,1.013,1.029,1.016,1.024,1.021,1.023];


cheetah10=[.223,.222,.223,.223,.225,.223,.225,.226,.227];
cheetah15=[.502,.500,.503,.504,.507,.504,.502,.504,.501];
cheetah20=[.897,.901,.902,.910,.908,.905,.905,.904,.902];
cheetah25=[1.414,1.410,1.408,1.411,1.416,1.409,1.413,1.412,1.403];
cheetah30=[1.949,1.943,1.937,1.941,1.946,1.934,1.932,1.941,1.938];

cheetahSpray10=[.230,.229,.229,.230,.232,.233,.234,.234,.232];
cheetahSpray15=[.537,.539,.541,.538,.538,.534,.530,.539,.536];
cheetahSpray20=[.951,.943,.945,.942,.940,.946,.949,.943,.944];
cheetahSpray25=[1.477,1.483,1.482,1.475,1.478,1.474,1.477,1.474,1.472];
cheetahSpray30=[2.070,2.048,2.043,2.041,2.044,2.037,2.046,2.042,2.043];



%the rig parameters
scaleToPivot=0.16;
pivotToStart=0.25;
startToEnd=0.195;
%%%%%%%%%%%%%%%%%%%%%Force in Newtons on the middle cylinder 
%convert the scale reading to the actual force on the cylinder, through the
%moment arm
FcylinderBlack10=(blackFur10(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderBlack15=(blackFur15(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderBlack20=(blackFur20(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderBlack25=(blackFur25(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderBlack30=(blackFur30(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);

FcylinderBlackSpray10=(blackSpray10(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderBlackSpray15=(blackSpray15(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderBlackSpray20=(blackSpray20(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderBlackSpray25=(blackSpray25(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderBlackSpray30=(blackSpray30(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);

FcylinderPipe10=(pipe10(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderPipe15=(pipe15(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderPipe20=(pipe20(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderPipe25=(pipe25(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcylinderPipe30=(pipe30(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);

Fcheetah10=(cheetah10(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
Fcheetah15=(cheetah15(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
Fcheetah20=(cheetah20(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
Fcheetah25=(cheetah25(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
Fcheetah30=(cheetah30(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);

FcheetahSpray10=(cheetahSpray10(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcheetahSpray15=(cheetahSpray15(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcheetahSpray20=(cheetahSpray20(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcheetahSpray25=(cheetahSpray25(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);
FcheetahSpray30=(cheetahSpray30(:)*scaleToPivot*9.81)/(pivotToStart+startToEnd/2);

%parameters for the pipe
diameterPipe=0.05;
diameterFur=0.052;
diameterSpray=0.052;
diameterCheetah=0.0835;%0.057;%0.0835 % 0.057% ;%%%%%%%%%effective diameter to make the result comprable 
diameterCheetahSpray=0.057;
%fur add 1cm each side 

%parameters based on the temerature for the day
% the 1 is for the first experiment (gorilla fur), the 2 is for the second experiment (cheetah fur)
abs2Pressure=101.016*1000;%1,010,16
absPressure=102.811*1000;%how it changes with altitude 600 to 800m 
Rspecific=287.058;
tempKelvin=290.35; 
temp2Kelvin=292.15; 
rho=absPressure/(Rspecific*tempKelvin);
rho2=abs2Pressure/(Rspecific*temp2Kelvin);
dynamicViscosity1=1.8054*10^(-5);
dynamicViscosity2=1.815*10^(-5);

%areas of the cylinders 
areaPipe=diameterPipe*startToEnd;
areaFur=diameterFur*startToEnd;
areaSpray=diameterSpray*startToEnd;
areaCheetah=diameterCheetah*startToEnd;
areaCheetahSpray=diameterCheetahSpray*startToEnd;


%calculate the drag coefficients 
CdBlackFur10=FcylinderBlack10(:)*2/(10*10*rho*areaFur);
CdBlackFur15=FcylinderBlack15(:)*2/(15*15*rho*areaFur);
CdBlackFur20=FcylinderBlack20(:)*2/(20*20*rho*areaFur);
CdBlackFur25=FcylinderBlack25(:)*2/(25*25*rho*areaFur);
CdBlackFur30=FcylinderBlack30(:)*2/(30*30*rho*areaFur);

CdBlackSpray10=FcylinderBlackSpray10(:)*2/(10*10*rho*areaSpray);
CdBlackSpray15=FcylinderBlackSpray15(:)*2/(15*15*rho*areaSpray);
CdBlackSpray20=FcylinderBlackSpray20(:)*2/(20*20*rho*areaSpray);
CdBlackSpray25=FcylinderBlackSpray25(:)*2/(25*25*rho*areaSpray);
CdBlackSpray30=FcylinderBlackSpray30(:)*2/(30*30*rho*areaSpray);

CdPipe10=FcylinderPipe10(:)*2/(10*10*rho*areaPipe);
CdPipe15=FcylinderPipe15(:)*2/(15*15*rho*areaPipe);
CdPipe20=FcylinderPipe20(:)*2/(20*20*rho*areaPipe);
CdPipe25=FcylinderPipe25(:)*2/(25*25*rho*areaPipe);
CdPipe30=FcylinderPipe30(:)*2/(30*30*rho*areaPipe);

CdCheetah10=Fcheetah10(:)*2/(10*10*rho2*areaCheetah);
CdCheetah15=Fcheetah15(:)*2/(15*15*rho2*areaCheetah);
CdCheetah20=Fcheetah20(:)*2/(20*20*rho2*areaCheetah);
CdCheetah25=Fcheetah25(:)*2/(25*25*rho2*areaCheetah);
CdCheetah30=Fcheetah30(:)*2/(30*30*rho2*areaCheetah);

CdCheetahSpray10=FcheetahSpray10(:)*2/(10*10*rho2*areaCheetahSpray);
CdCheetahSpray15=FcheetahSpray15(:)*2/(15*15*rho2*areaCheetahSpray);
CdCheetahSpray20=FcheetahSpray20(:)*2/(20*20*rho2*areaCheetahSpray);
CdCheetahSpray25=FcheetahSpray25(:)*2/(25*25*rho2*areaCheetahSpray);
CdCheetahSpray30=FcheetahSpray30(:)*2/(30*30*rho2*areaCheetahSpray);


%combine the vectors to generate a force vector(vector of all the force
%readings- makes plotting easier)
FfurVec(1,:)=[FcylinderBlack10(:)',FcylinderBlack15(:)',FcylinderBlack20(:)',FcylinderBlack25(:)',FcylinderBlack30(:)'];
FsprayVec(1,:)=[FcylinderBlackSpray10(:)',FcylinderBlackSpray15(:)',FcylinderBlackSpray20(:)',FcylinderBlackSpray25(:)',FcylinderBlackSpray30(:)'];
FpipeVec(1,:)=[FcylinderPipe10(:)',FcylinderPipe15(:)',FcylinderPipe20(:)',FcylinderPipe25(:)',FcylinderPipe30(:)'];
FcheetahVec(1,:)=[Fcheetah10(:)',Fcheetah15(:)',Fcheetah20(:)',Fcheetah25(:)',Fcheetah30(:)'];
FcheetahSprayVec(1,:)=[FcheetahSpray10(:)',FcheetahSpray15(:)',FcheetahSpray20(:)',FcheetahSpray25(:)',FcheetahSpray30(:)'];

speed=[10,10,10,10,10,10,10,10,10,15,15,15,15,15,15,15,15,15,20,20,20,20,20,20,20,20,20,25,25,25,25,25,25,25,25,25,30,30,30,30,30,30,30,30,30];

%comine thevectors to generate a Cd vector- makes plotting easier
CdFurVec(1,:)=[CdBlackFur10(:)',CdBlackFur15(:)',CdBlackFur20(:)',CdBlackFur25(:)',CdBlackFur30(:)'];
CdSprayVec(1,:)=[CdBlackSpray10(:)',CdBlackSpray15(:)',CdBlackSpray20(:)',CdBlackSpray25(:)',CdBlackSpray30(:)'];
CdPipeVec(1,:)=[CdPipe10(:)',CdPipe15(:)',CdPipe20(:)',CdPipe25(:)',CdPipe30(:)'];
CdCheetahVec(1,:)=[CdCheetah10(:)',CdCheetah15(:)',CdCheetah20(:)',CdCheetah25(:)',CdCheetah30(:)'];
CdCheetahSprayVec(1,:)=[CdCheetahSpray10(:)',CdCheetahSpray15(:)',CdCheetahSpray20(:)',CdCheetahSpray25(:)',CdCheetahSpray30(:)'];

%scatter plots
figure(1)
plot(speed(:),CdFurVec(:),'c+');%*3333 loglog
hold on
plot(speed(:),CdSprayVec(:),'k+');
plot(speed(:),CdPipeVec(:),'r+');
plot(speed(:),CdCheetahVec(:),'m+');
plot(speed(:),CdCheetahSprayVec(:),'b+');
xlabel('speed (m/s)');
ylabel('Cd');
title('Cd vs speed for all (vector) SCATTER');
legend('fur','spray','pipe','cheetah','cheetah spray');

figure(2)
plot(speed(:),FfurVec(:),'c+');
hold on
plot(speed(:),FsprayVec(:),'k+');
plot(speed(:),FpipeVec(:),'r+');
plot(speed(:),FcheetahVec(:),'m+');
plot(speed(:),FcheetahSprayVec(:),'b+');
xlabel('speed (m/s)');
ylabel('force');
title('force vs speed for all (vector) SCATTER');
legend('fur','spray','pipe','cheetah','cheetah spray');




%the error bar plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fvel1=[10,15,20,25,30];
y1=[mean(CdFurVec(1:9)),mean(CdFurVec(10:18)),mean(CdFurVec(19:27)),mean(CdFurVec(28:36)),mean(CdFurVec(37:45))];
y2=[mean(CdSprayVec(1:9)),mean(CdSprayVec(10:18)),mean(CdSprayVec(19:27)),mean(CdSprayVec(28:36)),mean(CdSprayVec(37:45))];
y3=[mean(CdPipeVec(1:9)),mean(CdPipeVec(10:18)),mean(CdPipeVec(19:27)),mean(CdPipeVec(28:36)),mean(CdPipeVec(37:45))];
y4=[mean(CdCheetahVec(1:9)),mean(CdCheetahVec(10:18)),mean(CdCheetahVec(19:27)),mean(CdCheetahVec(28:36)),mean(CdCheetahVec(37:45))];
y5=[mean(CdCheetahSprayVec(1:9)),mean(CdCheetahSprayVec(10:18)),mean(CdCheetahSprayVec(19:27)),mean(CdCheetahSprayVec(28:36)),mean(CdCheetahSprayVec(37:45))];

u1=[max(CdFurVec(1:9)),max(CdFurVec(10:18)),max(CdFurVec(19:27)),max(CdFurVec(28:36)),max(CdFurVec(37:45))]-y1;
u2=[max(CdSprayVec(1:9)),max(CdSprayVec(10:18)),max(CdSprayVec(19:27)),max(CdSprayVec(28:36)),max(CdSprayVec(37:45))]-y2;
u3=[max(CdPipeVec(1:9)),max(CdPipeVec(10:18)),max(CdPipeVec(19:27)),max(CdPipeVec(28:36)),max(CdPipeVec(37:45))]-y3;
u4=[max(CdCheetahVec(1:9)),max(CdCheetahVec(10:18)),max(CdCheetahVec(19:27)),max(CdCheetahVec(28:36)),max(CdCheetahVec(37:45))]-y4;
u5=[max(CdCheetahSprayVec(1:9)),max(CdCheetahSprayVec(10:18)),max(CdCheetahSprayVec(19:27)),max(CdCheetahSprayVec(28:36)),max(CdCheetahSprayVec(37:45))]-y5;

l1=y1-[min(CdFurVec(1:9)),min(CdFurVec(10:18)),min(CdFurVec(19:27)),min(CdFurVec(28:36)),min(CdFurVec(37:45))];
l2=y2-[min(CdSprayVec(1:9)),min(CdSprayVec(10:18)),min(CdSprayVec(19:27)),min(CdSprayVec(28:36)),min(CdSprayVec(37:45))];
l3=y3-[min(CdPipeVec(1:9)),min(CdPipeVec(10:18)),min(CdPipeVec(19:27)),min(CdPipeVec(28:36)),min(CdPipeVec(37:45))];
l4=y4-[min(CdCheetahVec(1:9)),min(CdCheetahVec(10:18)),min(CdCheetahVec(19:27)),min(CdCheetahVec(28:36)),min(CdCheetahVec(37:45))];
l5=y5-[min(CdCheetahSprayVec(1:9)),min(CdCheetahSprayVec(10:18)),min(CdCheetahSprayVec(19:27)),min(CdCheetahSprayVec(28:36)),min(CdCheetahSprayVec(37:45))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scatter plots
figure(3)
%errorbar(Fvel1(:),y1(:),l1,u1,'c+');
hold on
%errorbar(Fvel1(:),y2(:),l2,u2,'k+');
%errorbar(Fvel1(:),y3(:),l3,u3,'r+');
%errorbar(Fvel1(:),y4(:),l4,u4,'m+');
%errorbar(Fvel1(:),y5(:),l5,u5,'b+');
% for i=1:length(y4)
%     y4(i)=y4(i)/
%plot(Fvel1(:),y1(:),'c');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%changed the diameter
%plot(Fvel1(:),y2(:),'k');
plot(Fvel1(:),y3(:),'Color','r','LineWidth',2);
plot(Fvel1(:),y4(:),'Color','b','LineWidth',2);
%plot(Fvel1(:),y5(:),'b');
D=diameterCheetah;
myData=y4;
xlabel('Airspeed (m/s)','FontSize',16);
ylabel('Drag Coefficient (Cd)','FontSize',14);
title('Drag Coefficient vs Airspeed','FontSize',14);
legend('PVC Pipe','Cheetah Fur');
box on
grid on

%plot for reynolds numbers 
%Cd vs reynolds numbers 
figure(4)
r1=[10,15,20,25,30]*rho*diameterFur/dynamicViscosity1;%fur
r2=[10,15,20,25,30]*rho*diameterPipe/dynamicViscosity1;%wood
r3=[10,15,20,25,30]*rho2*diameterCheetah/dynamicViscosity2;%cheetah
semilogx(r1,y1,'c');
hold on
semilogx(r1,y2,'k');
semilogx(r2,y3,'r');
semilogx(r3,y4,'m');
semilogx(r3,y5,'b');
xlabel('Reynolds number');
ylabel('Cd');
title('Cd vs reynolds number');
legend('fur','spray','pipe','cheetah','cheetah spray');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fvel1=[10,15,20,25,30];
y1=[mean(FfurVec(1:9)),mean(FfurVec(10:18)),mean(FfurVec(19:27)),mean(FfurVec(28:36)),mean(FfurVec(37:45))];
y2=[mean(FsprayVec(1:9)),mean(FsprayVec(10:18)),mean(FsprayVec(19:27)),mean(FsprayVec(28:36)),mean(FsprayVec(37:45))];
y3=[mean(FpipeVec(1:9)),mean(FpipeVec(10:18)),mean(FpipeVec(19:27)),mean(FpipeVec(28:36)),mean(FpipeVec(37:45))];
y4=[mean(FcheetahVec(1:9)),mean(FcheetahVec(10:18)),mean(FcheetahVec(19:27)),mean(FcheetahVec(28:36)),mean(FcheetahVec(37:45))];
y5=[mean(FcheetahSprayVec(1:9)),mean(FcheetahSprayVec(10:18)),mean(FcheetahSprayVec(19:27)),mean(FcheetahSprayVec(28:36)),mean(FcheetahSprayVec(37:45))];

u1=[max(FfurVec(1:9)),max(FfurVec(10:18)),max(FfurVec(19:27)),max(FfurVec(28:36)),max(FfurVec(37:45))]-y1;
u2=[max(FsprayVec(1:9)),max(FsprayVec(10:18)),max(FsprayVec(19:27)),max(FsprayVec(28:36)),max(FsprayVec(37:45))]-y2;
u3=[max(FpipeVec(1:9)),max(FpipeVec(10:18)),max(FpipeVec(19:27)),max(FpipeVec(28:36)),max(FpipeVec(37:45))]-y3;
u4=[max(FcheetahVec(1:9)),max(FcheetahVec(10:18)),max(FcheetahVec(19:27)),max(FcheetahVec(28:36)),max(FcheetahVec(37:45))]-y4;
u5=[max(FcheetahSprayVec(1:9)),max(FcheetahSprayVec(10:18)),max(FcheetahSprayVec(19:27)),max(FcheetahSprayVec(28:36)),max(FcheetahSprayVec(37:45))]-y5;

l1=y1-[min(FfurVec(1:9)),min(FfurVec(10:18)),min(FfurVec(19:27)),min(FfurVec(28:36)),min(FfurVec(37:45))];
l2=y2-[min(FsprayVec(1:9)),min(FsprayVec(10:18)),min(FsprayVec(19:27)),min(FsprayVec(28:36)),min(FsprayVec(37:45))];
l3=y3-[min(FpipeVec(1:9)),min(FpipeVec(10:18)),min(FpipeVec(19:27)),min(FpipeVec(28:36)),min(FpipeVec(37:45))];
l4=y4-[min(FcheetahVec(1:9)),min(FcheetahVec(10:18)),min(FcheetahVec(19:27)),min(FcheetahVec(28:36)),min(FcheetahVec(37:45))];
l5=y5-[min(FcheetahSprayVec(1:9)),min(FcheetahSprayVec(10:18)),min(FcheetahSprayVec(19:27)),min(FcheetahSprayVec(28:36)),min(FcheetahSprayVec(37:45))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure(5)
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
ylabel('force');
title('force vs speed for all (vector) SCATTER');
legend('fur','spray','pipe','cheetah','cheetah spray');





%Reynolds numbers 
%1.8-dynaic viscosity -temp dep 1.983 x 10-5
%Re=density*vel*diameter/dynamicViscosity;





