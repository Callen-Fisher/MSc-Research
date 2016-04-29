function [x,f] = runtracklsq(x0,myfile)
% RUNTRACKLSQ demonstrates using LSQNONLIN with Simulink.

initCond=x0;

global mydata;
mydata=load(myfile);

options = optimoptions(@fmincon,'Algorithm','interior-point','display','iter','MaxFunEvals',15000,'MaxIter',7000,...
    'FinDiffRelStep',1e-4, 'TolFun', 1e-4, 'TolX', 1e-4,'diffminchange',1e-4,'UseParallel','always');
warning('off');
problem = createOptimProblem('fmincon','x0',initCond,'objective',@tracklsq,...
    'lb',[0.0000000000001,0.0000000000001]','ub',[1,1]','options',options);
gs = GlobalSearch('Display','iter','NumTrialPoints',1500);
[x,f]=run(gs,problem);

    function F = tracklsq(vals)
      
      mydata;
        sampleTime=mydata.sampleTime;
        states=mydata.states;
        val1 =mydata.val1;
        val2=mydata.val2;
        g1=mydata.g1;
        g2=mydata.g2;
        covP=mydata.covP;
        camAngles=mydata.camAngles;
        a2=mydata.a2;
        a1=mydata.a1;
        I=mydata.I;
        R=mydata.R;
        ma1=mydata.ma1;
        ma2=mydata.ma2;
        scalFac=1;
        l1=mydata.l1;
        l2=mydata.l2;
        camData=mydata.camData;
      val3 = vals(1)*scalFac;
      val4 = vals(2)*scalFac;

      
      [F,temp,temp2,temp3]=fullEKF(val1,val2,val3,val4,states,I,R,a1,a2,sampleTime,g1,g2,covP,camAngles,ma1,ma2,l1,l2,camData);
      
    end
end