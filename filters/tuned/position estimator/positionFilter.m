function [ costFunction,storedStates,storedPosition ] = positionFilter(acc,mag,gyro,g,magVecX,magVecY,magVecZ,R1,R2,I,covP,val1,val2,val3,val4,val5,val6,truePosition,states,sampleTime,gps)

Q=diag([val1,val2,val3,val4,val5,val6,0,0,0]);
startTime=1;
counter=0;
counterGPS=1;
for i=startTime:1:length(acc.signals.values(:,1)) 
    
    Rr=[[1 0 0];[0 cos(states(1)) -sin(states(1))];[0 sin(states(1)) cos(states(1))]];
    Rp=[[cos(states(2)) 0 sin(states(2))];[0 1 0];[-sin(states(2)) 0 cos(states(2))]];
    Ry=[[cos(states(3)) -sin(states(3)) 0];[sin(states(3)) cos(states(3)) 0];[0 0 1]];

    R=Ry*Rp*Rr;
    temp=R*[acc.signals.values(i,1);acc.signals.values(i,2);acc.signals.values(i,3)]-[0;0;g];
    ax=temp(1);
    ay=temp(2);
    az=temp(3);
    
    R_gyro=[[1 0 sin(states(2))];[0 cos(states(1)) -sin(states(1))*cos(states(2))];[0 sin(states(1)) cos(states(1))*cos(states(2))]];
    temp=R_gyro*[gyro.signals.values(i,3);gyro.signals.values(i,1);gyro.signals.values(i,2)];
    gx=temp(1);
    gy=temp(2);
    gz=temp(3);

    predictEq=predictFunc(ax,ay,az,0,gx,gy,gz,states(3),states(1),states(2),states(4),states(5),states(6));
   
    Fmatrix=Ffunc(ax,ay,az,gy,gz,states(3),states(1),states(2));
    Fmatrix=I+sampleTime*Fmatrix;
    
    states = states+sampleTime*predictEq';
    
    covP = (Fmatrix)*(covP*(Fmatrix')) + Q;
    counter=0;
    if(counter==10)
        counter=0;
        updateEq2=updateFunc2(acc.signals.values(i,1),acc.signals.values(i,2),acc.signals.values(i,3),gyro.signals.values(i,3),gyro.signals.values(i,1),gyro.signals.values(i,2),g,magVecX,magVecY,magVecZ,states(3),states(7),states(8),states(9),states(1),states(2));
        Hmatrix2=Hfunc2(acc.signals.values(i,1),acc.signals.values(i,2),acc.signals.values(i,3),gyro.signals.values(i,1),gyro.signals.values(i,2),g,magVecX,magVecY,magVecZ,states(3),states(1),states(2));
        K2=(covP*Hmatrix2')/(Hmatrix2*covP*Hmatrix2'+R2);
        z2=[acc.signals.values(i,1) acc.signals.values(i,2) acc.signals.values(i,3) gyro.signals.values(i,1) gyro.signals.values(i,2) gyro.signals.values(i,3) mag.signals.values(i,1) mag.signals.values(i,2) mag.signals.values(i,3),gps.signals.values(counterGPS,1),gps.signals.values(counterGPS,2),gps.signals.values(counterGPS,3)]'; 
        states=states+K2*(z2-updateEq2');
        covP=(I-K2*Hmatrix2)*covP;
        counterGPS=counterGPS+1;
    else
        updateEq1=updateFunc1(acc.signals.values(i,1),acc.signals.values(i,2),acc.signals.values(i,3),gyro.signals.values(i,3),gyro.signals.values(i,1),gyro.signals.values(i,2),g,magVecX,magVecY,magVecZ,states(3),states(1),states(2));
        Hmatrix1=Hfunc1(acc.signals.values(i,1),acc.signals.values(i,2),acc.signals.values(i,3),gyro.signals.values(i,1),gyro.signals.values(i,2),g,magVecX,magVecY,magVecZ,states(3),states(1),states(2));
        K1=(covP*Hmatrix1')/(Hmatrix1*covP*Hmatrix1'+R1);
        z1=[acc.signals.values(i,1) acc.signals.values(i,2) acc.signals.values(i,3) gyro.signals.values(i,1) gyro.signals.values(i,2) gyro.signals.values(i,3) mag.signals.values(i,1) mag.signals.values(i,2) mag.signals.values(i,3)]'; 
        states=states+K1*(z1-updateEq1');
        covP=(I-K1*Hmatrix1)*covP;
    end
    storedStates(i,:)=states(:)';
    storedPosition(i,:)=[states(7),states(8),states(9)];
    
    counter=counter+1;
    
    
end

%%%%%%%%%%%%%%%cam data stuff 
costFunction=sqrt((sum((storedPosition(:,1)-truePosition.signals.values(:,1)).^2)+...
                   sum((storedPosition(:,2)-truePosition.signals.values(:,2)).^2)+...
                   sum((storedPosition(:,3)-truePosition.signals.values(:,3)).^2))/(length(acc.signals.values(:,1))));
end

