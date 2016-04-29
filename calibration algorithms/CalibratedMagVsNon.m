clc;clear

[ U,c ] = calibrateMagFunction( 'mag.dat' );

load('rawMag')

figure(1)

plot3(magS1(:,1),magS1(:,2),magS1(:,3),'LineWidth',2);
xlabel('X axis','FontSize',10)
ylabel('Y axis','FontSize',10)
zlabel('Z axis','FontSize',10)
title('Raw, Uncalibrated Magnetometer Data','FontSize',14)




figure(4)
plot(magS1(:,1)+0.1,magS1(:,3)*0.6,'LineWidth',2);
xlabel('X axis','FontSize',14)
ylabel('Z axis','FontSize',14)
axis([-0.5 0.5 -0.5 0.5])
title('Raw, Uncalibrated Magnetometer Data','FontSize',18)






%calibrate the data
U=    [[3.5589    0.0788   -0.0108]
         [0    4.4166    0.0803]
         [0         0    3.1110]];
for i=1:1:length(magS1(:,1))
    w(i,:) = U*(magS1(i,:)'-[0.0003 0.0058 -0.1819]');
end

figure(2)
plot3(w(:,1),w(:,2),w(:,3),'LineWidth',2);
xlabel('X axis','FontSize',10)
ylabel('Y axis','FontSize',10)
zlabel('Z axis','FontSize',10)
title('Calibrated Magnetometer Data','FontSize',14)




% figure(3)
% plot3(magS1(:,1),magS1(:,2),magS1(:,3),'LineWidth',2);
% hold on
% plot3(w(:,1),w(:,2),w(:,3),'r','LineWidth',2);
% 
% legend('uncalibrated','calibrated')
% title('Calibrated Versus Uncalibrated Magnetometer Data','FontSize',14)



figure(5)
plot(w(:,1),w(:,3),'LineWidth',2);
xlabel('X axis','FontSize',14)
ylabel('Z axis','FontSize',14)
axis([-1 1 -1 1])
title('Calibrated Magnetometer Data','FontSize',18)