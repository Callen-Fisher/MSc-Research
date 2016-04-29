function [x,f] = runtracklsq(x0,myfile)
% RUNTRACKLSQ demonstrates using LSQNONLIN with Simulink.

initCond=x0;

global mydata;
mydata=load(myfile);

options = optimoptions(@fmincon,'Algorithm','interior-point','display','iter','MaxFunEvals',15000,'MaxIter',7000,...
    'FinDiffRelStep',1e-4, 'TolFun', 1e-4, 'TolX', 1e-4,'diffminchange',1e-4,'UseParallel','always');
warning('off');
problem = createOptimProblem('fmincon','x0',initCond,'objective',@tracklsq,...
    'lb',[0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001,0.0000000000001]','ub',[1,1,1,1,1,1,1,1,1,1]','options',options);
gs = GlobalSearch('Display','iter','NumTrialPoints',1500);
[x,f]=run(gs,problem);


    function F = tracklsq(vals)
      
      mydata;
        sampleTime=mydata.sampleTime;
        states=mydata.states;
        g1=mydata.g1;
        g2=mydata.g2;
        g3=mydata.g3;
        g4=mydata.g4;
        covP=mydata.covP;
        camData=mydata.camData;
        a4=mydata.a4;
        a3=mydata.a3;
        a2=mydata.a2;
        a1=mydata.a1;
        I=mydata.I;
        R=mydata.R;
        l1=mydata.l1;
        l2=mydata.l2;
        l3=mydata.l3;
        l4=mydata.l4;
        
        ma1=mydata.ma1;
        ma2=mydata.ma2;
        ma3=mydata.ma3;
        ma4=mydata.ma4;
        scalFac=1;
        
      val3 = vals(1)*scalFac;
      val4 = vals(2)*scalFac;
      val5 =vals(3)*scalFac;
      val6=vals(4)*scalFac;
      val7=vals(5)*scalFac;
      val8=vals(6)*scalFac;
      val9=vals(7)*scalFac;
      val10=vals(8)*scalFac;
      val11=vals(9)*scalFac;
      val12=vals(10)*scalFac;
      
      [F,temp,temp2]=EKF(val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,states,covP,R,I,sampleTime,camData,a1,a2,a3,a4,g1,g2,g3,g4,ma1,ma2,ma3,ma4,l1,l2,l3,l4);
      
    end
end