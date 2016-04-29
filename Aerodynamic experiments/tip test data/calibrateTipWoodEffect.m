
%HACK

scaleToPivot=0.16*9.81;%%%%%%%%%%%%HACK


pivotToStart=0.12;
lengthWood=0.035;
a=0.024-2*0.0008873;%width of smaller edge
b=0.024;%width of bigger edge 
x=(lengthWood/3)*(2*a+b)/(a+b);%distance to centre of area
area=(a+b)*lengthWood/2;
Rspecific=287.058;
absPressure=101.829*1000;
tempKelvin=291.45;
rho=absPressure/(Rspecific*tempKelvin);
cd10=1.0332;
cd15=1.0672;
cd20=1.0432;
cd25=0.8884;
cd30=0.8508;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd=cd10;
Fd=1/2*rho*area*cd*10^2;

%%%%%HACK
pivotToStart=0.12;
x=x;

blackFurTip10_0(:)=blackFurTip10_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;%first cos is for the length x and the second is for the area
blackFurTip10_15(:)=blackFurTip10_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackFurTip10_30(:)=blackFurTip10_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackFurTip10_45(:)=blackFurTip10_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackFurTip10_60(:)=blackFurTip10_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd15;
Fd=1/2*rho*area*cd*15^2;

blackFurTip15_0(:)=blackFurTip15_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
blackFurTip15_15(:)=blackFurTip15_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackFurTip15_30(:)=blackFurTip15_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackFurTip15_45(:)=blackFurTip15_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackFurTip15_60(:)=blackFurTip15_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd20;
Fd=1/2*rho*area*cd*20^2;

blackFurTip20_0(:)=blackFurTip20_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
blackFurTip20_15(:)=blackFurTip20_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackFurTip20_30(:)=blackFurTip20_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackFurTip20_45(:)=blackFurTip20_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackFurTip20_60(:)=blackFurTip20_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd25;
Fd=1/2*rho*area*cd*25^2;

blackFurTip25_0(:)=blackFurTip25_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
blackFurTip25_15(:)=blackFurTip25_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackFurTip25_30(:)=blackFurTip25_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackFurTip25_45(:)=blackFurTip25_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackFurTip25_60(:)=blackFurTip25_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd30;
Fd=1/2*rho*area*cd*30^2;

blackFurTip30_0(:)=blackFurTip30_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
blackFurTip30_15(:)=blackFurTip30_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackFurTip30_30(:)=blackFurTip30_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackFurTip30_45(:)=blackFurTip30_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackFurTip30_60(:)=blackFurTip30_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd=cd10;
Fd=1/2*rho*area*cd*10^2;

blackSprayTip10_0(:)=blackSprayTip10_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
blackSprayTip10_15(:)=blackSprayTip10_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackSprayTip10_30(:)=blackSprayTip10_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackSprayTip10_45(:)=blackSprayTip10_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackSprayTip10_60(:)=blackSprayTip10_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd15;
Fd=1/2*rho*area*cd*15^2;

blackSprayTip15_0(:)=blackSprayTip15_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
blackSprayTip15_15(:)=blackSprayTip15_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackSprayTip15_30(:)=blackSprayTip15_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackSprayTip15_45(:)=blackSprayTip15_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackSprayTip15_60(:)=blackSprayTip15_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd20;
Fd=1/2*rho*area*cd*20^2;

blackSprayTip20_0(:)=blackSprayTip20_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
blackSprayTip20_15(:)=blackSprayTip20_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackSprayTip20_30(:)=blackSprayTip20_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackSprayTip20_45(:)=blackSprayTip20_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackSprayTip20_60(:)=blackSprayTip20_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd25;
Fd=1/2*rho*area*cd*25^2;

blackSprayTip25_0(:)=blackSprayTip25_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
blackSprayTip25_15(:)=blackSprayTip25_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackSprayTip25_30(:)=blackSprayTip25_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackSprayTip25_45(:)=blackSprayTip25_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackSprayTip25_60(:)=blackSprayTip25_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd30;
Fd=1/2*rho*area*cd*30^2;

blackSprayTip30_0(:)=blackSprayTip30_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
blackSprayTip30_15(:)=blackSprayTip30_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
blackSprayTip30_30(:)=blackSprayTip30_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
blackSprayTip30_45(:)=blackSprayTip30_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
blackSprayTip30_60(:)=blackSprayTip30_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lengthWood=0.10;
a=0.01893;%width of smaller edge
b=0.024;%width of bigger edge 
x=(lengthWood/3)*(2*a+b)/(a+b);%distance to centre of area
area=(a+b)*lengthWood/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%HACK
pivotToStart=0.12;
x=x;



cd=cd10;
Fd=1/2*rho*area*cd*10^2;

cheetahSprayTip10_0(:)=cheetahSprayTip10_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahSprayTip10_15(:)=cheetahSprayTip10_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahSprayTip10_30(:)=cheetahSprayTip10_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahSprayTip10_45(:)=cheetahSprayTip10_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahSprayTip10_60(:)=cheetahSprayTip10_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd15;
Fd=1/2*rho*area*cd*15^2;

cheetahSprayTip15_0(:)=cheetahSprayTip15_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahSprayTip15_15(:)=cheetahSprayTip15_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahSprayTip15_30(:)=cheetahSprayTip15_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahSprayTip15_45(:)=cheetahSprayTip15_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahSprayTip15_60(:)=cheetahSprayTip15_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd20;
Fd=1/2*rho*area*cd*20^2;

cheetahSprayTip20_0(:)=cheetahSprayTip20_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahSprayTip20_15(:)=cheetahSprayTip20_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahSprayTip20_30(:)=cheetahSprayTip20_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahSprayTip20_45(:)=cheetahSprayTip20_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahSprayTip20_60(:)=cheetahSprayTip20_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd25;
Fd=1/2*rho*area*cd*25^2;

cheetahSprayTip25_0(:)=cheetahSprayTip25_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahSprayTip25_15(:)=cheetahSprayTip25_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahSprayTip25_30(:)=cheetahSprayTip25_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahSprayTip25_45(:)=cheetahSprayTip25_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahSprayTip25_60(:)=cheetahSprayTip25_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd30;
Fd=1/2*rho*area*cd*30^2;

cheetahSprayTip30_0(:)=cheetahSprayTip30_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahSprayTip30_15(:)=cheetahSprayTip30_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahSprayTip30_30(:)=cheetahSprayTip30_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahSprayTip30_45(:)=cheetahSprayTip30_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahSprayTip30_60(:)=cheetahSprayTip30_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd=cd10;
Fd=1/2*rho*area*cd*10^2;

cheetahTip10_0(:)=cheetahTip10_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahTip10_15(:)=cheetahTip10_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahTip10_30(:)=cheetahTip10_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahTip10_45(:)=cheetahTip10_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahTip10_60(:)=cheetahTip10_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd15;
Fd=1/2*rho*area*cd*15^2;

cheetahTip15_0(:)=cheetahTip15_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahTip15_15(:)=cheetahTip15_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahTip15_30(:)=cheetahTip15_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahTip15_45(:)=cheetahTip15_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahTip15_60(:)=cheetahTip15_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd20;
Fd=1/2*rho*area*cd*20^2;

cheetahTip20_0(:)=cheetahTip20_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahTip20_15(:)=cheetahTip20_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahTip20_30(:)=cheetahTip20_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahTip20_45(:)=cheetahTip20_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahTip20_60(:)=cheetahTip20_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd25;
Fd=1/2*rho*area*cd*25^2;

cheetahTip25_0(:)=cheetahTip25_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahTip25_15(:)=cheetahTip25_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahTip25_30(:)=cheetahTip25_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahTip25_45(:)=cheetahTip25_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahTip25_60(:)=cheetahTip25_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;

cd=cd30;
Fd=1/2*rho*area*cd*30^2;

cheetahTip30_0(:)=cheetahTip30_0(:)-Fd*(pivotToStart+x*cos(0*pi/180))*cos(0*pi/180)/scaleToPivot;
cheetahTip30_15(:)=cheetahTip30_15(:)-Fd*(pivotToStart+x*cos(15*pi/180))*cos(15*pi/180)/scaleToPivot;
cheetahTip30_30(:)=cheetahTip30_30(:)-Fd*(pivotToStart+x*cos(30*pi/180))*cos(30*pi/180)/scaleToPivot;
cheetahTip30_45(:)=cheetahTip30_45(:)-Fd*(pivotToStart+x*cos(45*pi/180))*cos(45*pi/180)/scaleToPivot;
cheetahTip30_60(:)=cheetahTip30_60(:)-Fd*(pivotToStart+x*cos(60*pi/180))*cos(60*pi/180)/scaleToPivot;