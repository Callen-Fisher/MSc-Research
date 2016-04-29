%key:
%camera=red
%2Dekf=blue
%2Dfull=green
%3Dekf=cyan
%3Dfull=black
%2Dtriad=magenta
%3Dtriad=yellow

%%%%%%%%%%%%Positions 
figure(1)
subplot(2,3,1)
plot(camData.signals.values(:,1)','r');
hold on
plot(storedPositions2Dekf(:,1)','b');
plot(storedPositions2Dfull(:,1)','g');
plot(storedPositions3Dekf(:,1)','c');
plot(storedPositions3Dfull(:,1)','k');
plot(storedPositions2Dtriad(:,1)','m');
plot(storedPositions3Dtriad(:,1)','y');

subplot(2,3,2)
plot(camData.signals.values(:,2)','r');
hold on
plot(storedPositions2Dekf(:,2)','b');
plot(storedPositions2Dfull(:,2)','g');
plot(storedPositions3Dekf(:,2)','c');
plot(storedPositions3Dfull(:,2)','k');
plot(storedPositions2Dtriad(:,2)','m');
plot(storedPositions3Dtriad(:,2)','y');

subplot(2,3,3)
plot(camData.signals.values(:,3)','r');
hold on
plot(storedPositions2Dekf(:,3)','b');
plot(storedPositions2Dfull(:,3)','g');
plot(storedPositions3Dekf(:,3)','c');
plot(storedPositions3Dfull(:,3)','k');
plot(storedPositions2Dtriad(:,3)','m');
plot(storedPositions3Dtriad(:,3)','y');

subplot(2,3,4)
plot(camData.signals.values(:,4)','r');
hold on
plot(storedPositions2Dekf(:,4)','b');
plot(storedPositions2Dfull(:,4)','g');
plot(storedPositions3Dekf(:,4)','c');
plot(storedPositions3Dfull(:,4)','k');
plot(storedPositions2Dtriad(:,4)','m');
plot(storedPositions3Dtriad(:,4)','y');

subplot(2,3,5)
plot(camData.signals.values(:,5)','r');
hold on
plot(storedPositions2Dekf(:,5)','b');
plot(storedPositions2Dfull(:,5)','g');
plot(storedPositions3Dekf(:,5)','c');
plot(storedPositions3Dfull(:,5)','k');
plot(storedPositions2Dtriad(:,5)','m');
plot(storedPositions3Dtriad(:,5)','y');

subplot(2,3,6)
plot(camData.signals.values(:,6)','r');
hold on
plot(storedPositions2Dekf(:,6)','b');
plot(storedPositions2Dfull(:,6)','g');
plot(storedPositions3Dekf(:,6)','c');
plot(storedPositions3Dfull(:,6)','k');
plot(storedPositions2Dtriad(:,6)','m');
plot(storedPositions3Dtriad(:,6)','y');
%%%%%%%%%%%%Angles    
figure(2)
subplot(2,2,1)
plot(camAngles(:,1)','r');
hold on
plot(storedStates2Dekf(:,3)','b');
plot(storedStates2Dfull(:,5)','g');
plot(storedStates3Dekf(:,5)','c');
plot(storedStates3Dfull(:,9)','k');
plot(storedStates2Dtriad(:,1)','m');
plot(storedStates3Dtriad(:,1)','y');

subplot(2,2,2)
plot(camAngles(:,2)','r');
hold on
plot(storedStates3Dekf(:,6)','c');
plot(storedStates3Dfull(:,10)','k');
plot(storedStates3Dtriad(:,2)','y');

subplot(2,2,3)
plot(camAngles(:,3)','r');
hold on
plot(storedStates2Dekf(:,4)','b');
plot(storedStates2Dfull(:,6)','g');
plot(storedStates3Dekf(:,7)','c');
plot(storedStates3Dfull(:,11)','k');
plot(storedStates2Dtriad(:,2)','m');
plot(storedStates3Dtriad(:,3)','y');

subplot(2,2,4)
plot(camAngles(:,4)','r');
hold on
plot(storedStates3Dekf(:,8)','c');
plot(storedStates3Dfull(:,12)','k');
plot(storedStates3Dtriad(:,4)','y');

%%%%%%%%%%%%Rates     
figure(3)
subplot(2,2,1)
plot(storedStates2Dekf(:,1)','b');
hold on
plot(storedStates2Dfull(:,3)','g');
plot(storedStates3Dekf(:,1)','c');
plot(storedStates3Dfull(:,5)','k');

subplot(2,2,2)
plot(storedStates3Dekf(:,2)','c');
hold on
plot(storedStates3Dfull(:,6)','k');

subplot(2,2,3)
plot(storedStates2Dekf(:,2)','b');
hold on
plot(storedStates2Dfull(:,4)','g');
plot(storedStates3Dekf(:,3)','c');
plot(storedStates3Dfull(:,7)','k');

subplot(2,2,4)
plot(storedStates3Dekf(:,4)','c');
hold on
plot(storedStates3Dfull(:,8)','k');
%%%%%%%%%%%%torques   
figure(4)
subplot(2,2,1)
plot(storedStates2Dfull(:,1)','g');
hold on
plot(storedStates3Dfull(:,1)','k');

subplot(2,2,2)
plot(storedStates3Dfull(:,2)','k');

subplot(2,2,3)
plot(storedStates2Dfull(:,2)','g');
hold on
plot(storedStates3Dfull(:,3)','k');

subplot(2,2,4)
plot(storedStates3Dfull(:,4)','k');