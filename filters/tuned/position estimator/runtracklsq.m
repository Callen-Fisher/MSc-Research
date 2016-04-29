function [x,f] = runtracklsq(x0,myfile)
% RUNTRACKLSQ demonstrates using LSQNONLIN with Simulink.

initCond=x0;

global mydata;
mydata=load(myfile);

options = optimoptions(@fmincon,'Algorithm','interior-point','display','iter','MaxFunEvals',15000,'MaxIter',7000,...
    'FinDiffRelStep',1e-4, 'TolFun', 1e-4, 'TolX', 1e-4,'diffminchange',1e-4,'UseParallel','always');
warning('off');
problem = createOptimProblem('fmincon','x0',initCond,'objective',@tracklsq,...
    'lb',[0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001]','ub',[1,1,1,1,1,1]','options',options);
gs = GlobalSearch('Display','iter','NumTrialPoints',1500);
[x,f]=run(gs,problem);


    function F = tracklsq(vals)
      
        mydata;
        acc=mydata.acc;
        mag=mydata.mag;
        gyro=mydata.gyro;
        g=mydata.g;
        magVecX=mydata.magVecX;
        magVecY=mydata.magVecY;
        magVecZ=mydata.magVecZ;
        R1=mydata.R1;
        R2=mydata.R2;
        I=mydata.I;
        covP=mydata.covP;
        truePosition=mydata.truePosition;
        states=mydata.states;
        sampleTime=mydata.sampleTime;
        GPS_position=mydata.GPS_position;
        
        scalFac=1;
        
        val1 = vals(1)*scalFac;
        val2 = vals(2)*scalFac;
        val3 = vals(3)*scalFac;
        val4 = vals(4)*scalFac;
        val5 = vals(5)*scalFac;
        val6 = vals(6)*scalFac;
        
      
      [F,temp1,temp2]=positionFilter(acc,mag,gyro,g,magVecX,magVecY,magVecZ,R1,R2,I,covP,val1,val2,val3,val4,val5,val6,truePosition,states,sampleTime,GPS_position);
      
    end
end