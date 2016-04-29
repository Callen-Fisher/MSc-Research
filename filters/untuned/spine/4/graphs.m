%key:
%camera=red
%2Dekf=blue
%2Dfull=green
%3Dekf=cyan
%3Dfull=black

%%%%%%%%%%%%Positions -> 2D ekf+full, 3D ekf+full, camera data ->3*3
figure(1)
subplot(3,3,1)
plot(camData.signals.values(:,1+3)','r');
hold on
plot(storedPositions2Dekf(:,1+3)','b');
plot(storedPositions2Dfull(:,1+3)','g');
plot(storedPositions3Dekf(:,1+3)','c');
plot(storedPositions3Dfull(:,1+3)','k');

subplot(3,3,2)
plot(camData.signals.values(:,2+3)','r');
hold on
plot(storedPositions2Dekf(:,2+3)','b');
plot(storedPositions2Dfull(:,2+3)','g');
plot(storedPositions3Dekf(:,2+3)','c');
plot(storedPositions3Dfull(:,2+3)','k');

subplot(3,3,3)
plot(camData.signals.values(:,3+3)','r');
hold on
plot(storedPositions2Dekf(:,3+3)','b');
plot(storedPositions2Dfull(:,3+3)','g');
plot(storedPositions3Dekf(:,3+3)','c');
plot(storedPositions3Dfull(:,3+3)','k');

subplot(3,3,4)
plot(camData.signals.values(:,4+3)','r');
hold on
plot(storedPositions2Dekf(:,4+3)','b');
plot(storedPositions2Dfull(:,4+3)','g');
plot(storedPositions3Dekf(:,4+3)','c');
plot(storedPositions3Dfull(:,4+3)','k');

subplot(3,3,5)
plot(camData.signals.values(:,5+3)','r');
hold on
plot(storedPositions2Dekf(:,5+3)','b');
plot(storedPositions2Dfull(:,5+3)','g');
plot(storedPositions3Dekf(:,5+3)','c');
plot(storedPositions3Dfull(:,5+3)','k');

subplot(3,3,6)
plot(camData.signals.values(:,6+3)','r');
hold on
plot(storedPositions2Dekf(:,6+3)','b');
plot(storedPositions2Dfull(:,6+3)','g');
plot(storedPositions3Dekf(:,6+3)','c');
plot(storedPositions3Dfull(:,6+3)','k');

subplot(3,3,7)
plot(camData.signals.values(:,7+3)','r');
hold on
plot(storedPositions2Dekf(:,7+3)','b');
plot(storedPositions2Dfull(:,7+3)','g');
plot(storedPositions3Dekf(:,7+3)','c');
plot(storedPositions3Dfull(:,7+3)','k');

subplot(3,3,8)
plot(camData.signals.values(:,8+3)','r');
hold on
plot(storedPositions2Dekf(:,8+3)','b');
plot(storedPositions2Dfull(:,8+3)','g');
plot(storedPositions3Dekf(:,8+3)','c');
plot(storedPositions3Dfull(:,8+3)','k');

subplot(3,3,9)
plot(camData.signals.values(:,9+3)','r');
hold on
plot(storedPositions2Dekf(:,9+3)','b');
plot(storedPositions2Dfull(:,9+3)','g');
plot(storedPositions3Dekf(:,9+3)','c');
plot(storedPositions3Dfull(:,9+3)','k');
%%%%%%%%%%%%Angles    -> 2D ekf+full, camera data              ->1*4
figure(2)
subplot(4,1,1)
plot(camAngles(:,1)','r');
hold on
plot(storedStates2Dekf(:,5)','b');
plot(storedStates2Dfull(:,7)','g');

subplot(4,1,2)
plot(camAngles(:,2)','r');
hold on
plot(storedStates2Dekf(:,6)','b');
plot(storedStates2Dfull(:,8)','g');

subplot(4,1,3)
plot(camAngles(:,3)','r');
hold on
plot(storedStates2Dekf(:,7)','b');
plot(storedStates2Dfull(:,9)','g');

subplot(4,1,4)
plot(camAngles(:,4)','r');
hold on
plot(storedStates2Dekf(:,8)','b');
plot(storedStates2Dfull(:,10)','g');
%%%%%%%%%%%%Rates     -> 2D ekf+full, 3D ekf+full              ->3*4
figure(3)
subplot(4,3,1)
plot(storedStates3Dekf(:,1)','c');
hold on
plot(storedStates3Dfull(:,5)','k');
subplot(4,3,2)
plot(storedStates2Dekf(:,1)','b');
hold on
plot(storedStates2Dfull(:,3)','g');
plot(storedStates3Dekf(:,2)','c');
plot(storedStates3Dfull(:,6)','k');
subplot(4,3,3)
plot(storedStates3Dekf(:,3)','c');
hold on
plot(storedStates3Dfull(:,7)','k');

subplot(4,3,4)
plot(storedStates3Dekf(:,4)','c');
hold on
plot(storedStates3Dfull(:,8)','k');
subplot(4,3,5)
plot(storedStates2Dekf(:,2)','b');
hold on
plot(storedStates2Dfull(:,4)','g');
plot(storedStates3Dekf(:,5)','c');
plot(storedStates3Dfull(:,9)','k');
subplot(4,3,6)
plot(storedStates3Dekf(:,6)','c');
hold on
plot(storedStates3Dfull(:,10)','k');

subplot(4,3,8)
plot(storedStates2Dekf(:,3)','b');
hold on
plot(storedStates2Dfull(:,5)','g');
plot(storedStates3Dekf(:,7)','c');
plot(storedStates3Dfull(:,11)','k');
subplot(4,3,9)
plot(storedStates3Dekf(:,8)','c');
hold on
plot(storedStates3Dfull(:,12)','k');

subplot(4,3,11)
plot(storedStates2Dekf(:,4)','b');
hold on
plot(storedStates2Dfull(:,6)','g');
plot(storedStates3Dekf(:,9)','c');
plot(storedStates3Dfull(:,13)','k');
subplot(4,3,12)
plot(storedStates3Dekf(:,10)','c');
hold on
plot(storedStates3Dfull(:,14)','k');

%%%%%%%%%%%%torques   -> 2D full, 3D full                      ->2*2
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