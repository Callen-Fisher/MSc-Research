function [x,f] = runtracklsq(x0,myfile)
% RUNTRACKLSQ demonstrates using LSQNONLIN with Simulink.

initCond=x0;

global mydata;
mydata=load(myfile);

options = optimoptions(@fmincon,'Algorithm','interior-point','display','iter','MaxFunEvals',15000,'MaxIter',7000,...
    'FinDiffRelStep',1e-4, 'TolFun', 1e-4, 'TolX', 1e-4,'diffminchange',1e-4,'UseParallel','always');
warning('off');
problem = createOptimProblem('fmincon','x0',initCond,'objective',@tracklsq,...
    'lb',[0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001]','ub',[10,10,10,10,10,10,10,10]','options',options);
gs = GlobalSearch('Display','iter','NumTrialPoints',1500);
[x,f]=run(gs,problem);

    function F = tracklsq(vals)
      
      mydata;
        sampleTime=mydata.sampleTime;
        states=mydata.states;

        ma2=mydata.ma2;
        ma1=mydata.ma1;
        g1=mydata.g1;
        g2=mydata.g2;
        covP=mydata.covP;
        camAngles=mydata.camAngles;
        a2=mydata.a2;
        a1=mydata.a1;
        I=mydata.I;
        Q=mydata.Q;
        l1=mydata.l1;
        l2=mydata.l2;
        camData=mydata.camData;
        scalFac=1;
magVecX=mydata.magVecX;
magVecY=mydata.magVecY;
magVecZ=mydata.magVecZ;
g=mydata.g;
        
      val1 = vals(1)*scalFac;
      val2 = vals(2)*scalFac;
      val3 = vals(3)*scalFac;
      val4 = vals(4)*scalFac;
val5 = vals(5)*scalFac;
      val6 = vals(6)*scalFac;
      val7 = vals(7)*scalFac;
      val8 = vals(8)*scalFac;
      
      [temp,F,temp2,temp3]=TRIAD(g1,g2,a1,a2,ma1,ma2,magVecX,magVecY,magVecZ,g,Q,I,covP,val1,val2,val3,val4,sampleTime,camAngles,camData,l1,l2,states,val5,val6,val7,val8);
      
    end
end