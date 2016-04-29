clc
clear

load('positionSpine');

figure(1)
plot(storedPositions(:,10),storedPositions(:,12),'r','LineWidth',4);
legend('spine tip', 'Location', 'SouthEast');
xlabel('X axis');
ylabel('Y axis','FontSize',12,'FontWeight','bold','Color','r');
title('My graph')




%%%%%%%%%%%%SUBPLOT
figure % create new figure
subplot(2,2,1) % first subplot
plot(x,y1)
title('First subplot')
y2 = sin(2*x); % define y2

subplot(2,2,2) % second subplot
plot(x,y2)
title('Second subplot')

suptitle('Main title')



%%%%%%%%%%BACKGROUND COLOR

figure('Color',[0.8 0.8 0.8]);
whitebg
set(gcf,'color','white')