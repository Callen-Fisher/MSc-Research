clc;
clear
load('data');
textBook=[%[10000 1]
              %[19000 1.0316]
              %[28000 1.0478]
              [37000 1.0478]
              [46000 1.0642]
              [55000 1.0809]
              [64000 1.0809]
              [73000 1.0642]
              [82000 1.0157]
              [91000 1]
              [100000 0.9544]];
              %[109000 0.8169]
              %[118000 0.6078]
              %[127000 0.3065]
              %[136000 0.2177]
              %[145000 0.2246]
              %[154000 0.2465]];
         
plot(textBook(:,1),textBook(:,2),'LineWidth',2)
hold on
abs2Pressure=101.016*1000;%1,010,16
absPressure=102.811*1000;%how it changes with altitude 600 to 800m 
Rspecific=287.058;
tempKelvin=290.35; 
temp2Kelvin=292.15; 
rho=absPressure/(Rspecific*tempKelvin);
rho2=abs2Pressure/(Rspecific*temp2Kelvin);
dynamicViscosity1=1.8054*10^(-5);
dynamicViscosity2=1.815*10^(-5);
diameterPipe=0.052;
r2=[10,15,20,25,30]*rho2*diameterPipe/dynamicViscosity2;
plot(r2,y3,'r','LineWidth',2)
legend('Textbook data','Smooth cylinder test')
title('Known Textbook Data Compared to Smooth Cylinder Test','FontSize',12)
xlabel('Reynolds Number (Re)','FontSize',10)
ylabel('Drag Coefficient (Cd)','FontSize',10)