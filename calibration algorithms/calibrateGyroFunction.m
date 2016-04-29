function [ ] = calibrateGyroFunction( gyroFileName )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%the file format:
%GyroX1,GyroY1,GyroZ1,AccX1,AccY1,AccZ1,MagX1,MagY1,MagZ1,Temp1,Temp2,rollCommand,pitchCommand,yawCommand,rollEncoder,pitchEncoder,yawEncoder,yawVel

%the file names:
%Gyro.dat

%[Gx1 Gy1 Gz1]=[Gx Gy Gz]-[Bx By Bz]
%Y=X-B
%Y=calibrated data
%X=raw sensor data
%B=the sensor bias

[ gyroS1, magS1, accS1, tempS1, tempADCS1, yawEncoder,yawVelocity] = calibrationReadDataFunction( gyroFileName );

figure1 = figure;
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'all');
plot1=plot(tempADCS1(:),gyroS1(:,1),'Parent',axes1,'DisplayName','data1');
title('Gyro X, deg/s versus temp');
xlabel('temp');
ylabel('rad/s');
xdata1 = get(plot1, 'xdata');
ydata1 = get(plot1, 'ydata');
xdata1 = xdata1(:);
ydata1 = ydata1(:);
axesLimits1 = xlim(axes1);
xplot1 = linspace(axesLimits1(1), axesLimits1(2));
coeffs1 = cell(1,1);
fitResults1 = polyfit(xdata1, ydata1, 1);
yplot1 = polyval(fitResults1, xplot1);
fittypesArray1(1) = 2;
coeffs1{1} = fitResults1;
fitLine1 = plot(xplot1,yplot1,'DisplayName','   linear','Parent',axes1,'Tag','linear','Color',[1 0 0]);
setLineOrder(axes1, fitLine1, plot1);
showEquations(fittypesArray1, coeffs1, 2, axes1);
legend(axes1,'show');


figure2 = figure;
axes2 = axes('Parent',figure2);
box(axes2,'on');
hold(axes2,'all');
plot2=plot(tempADCS1(:),gyroS1(:,2),'Parent',axes2,'DisplayName','data2');
title('Gyro Y, deg/s versus temp');
xlabel('temp');
ylabel('rad/s');
xdata2 = get(plot2, 'xdata');
ydata2 = get(plot2, 'ydata');
xdata2 = xdata2(:);
ydata2 = ydata2(:);
axesLimits2 = xlim(axes2);
xplot2 = linspace(axesLimits2(1), axesLimits2(2));
coeffs2 = cell(1,1);
fitResults2 = polyfit(xdata2, ydata2, 1);
yplot2 = polyval(fitResults2, xplot2);
fittypesArray2(1) = 2;
coeffs2{1} = fitResults2;
fitLine2 = plot(xplot2,yplot2,'DisplayName','   linear','Parent',axes2,'Tag','linear','Color',[1 0 0]);
setLineOrder(axes2, fitLine2, plot2);
showEquations(fittypesArray2, coeffs2, 2, axes2);
legend(axes2,'show');


figure3 = figure;
axes3 = axes('Parent',figure3);
box(axes3,'on');
hold(axes3,'all');
plot3=plot(tempADCS1(:),gyroS1(:,3),'Parent',axes3,'DisplayName','data3');
title('Gyro Z, deg/s versus temp');
xlabel('temp');
ylabel('rad/s');
xdata3 = get(plot3, 'xdata');
ydata3 = get(plot3, 'ydata');
xdata3 = xdata3(:);
ydata3 = ydata3(:);
axesLimits3 = xlim(axes3);
xplot3 = linspace(axesLimits3(1), axesLimits3(2));
coeffs3 = cell(1,1);
fitResults3 = polyfit(xdata3, ydata3, 1);
yplot3 = polyval(fitResults3, xplot3);
fittypesArray3(1) = 2;
coeffs3{1} = fitResults3;
fitLine3 = plot(xplot3,yplot3,'DisplayName','   linear','Parent',axes3,'Tag','linear','Color',[1 0 0]);
setLineOrder(axes3, fitLine3, plot3);
showEquations(fittypesArray3, coeffs3, 2, axes3);
legend(axes3,'show');

end

