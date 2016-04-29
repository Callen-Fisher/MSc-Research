%code to calculate the force on the wood tail through the moment arm 
%parameters
scaleToPivot=0.16;
pivotToStart=0.12;
lengthTail=0.355;
a=0.006;%width of smaller edge
b=0.024;%width of bigger edge 
x=(lengthTail/3)*(2*a+b)/(a+b);%distance to centre of area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0=sqrt(pivotToStart^2+x^2-2*pivotToStart*x2*cos((180-0)*pi/180));%moment arm of actual force and not scale force
theta=acos((pivotToStartWood^2+c0^2-x2^2)/(2*pivotToStartWood*c0));
c0=c0*cos(theta);
c15=sqrt(pivotToStart^2+x^2-2*pivotToStart*x2*cos((180-15)*pi/180));
theta=acos((pivotToStartWood^2+c15^2-x2^2)/(2*pivotToStartWood*c15));
c15=c15*cos(theta);
c30=sqrt(pivotToStart^2+x^2-2*pivotToStart*x2*cos((180-30)*pi/180));
theta=acos((pivotToStartWood^2+c30^2-x2^2)/(2*pivotToStartWood*c30));
c30=c30*cos(theta);
c45=sqrt(pivotToStart^2+x^2-2*pivotToStart*x2*cos((180-45)*pi/180));
theta=acos((pivotToStartWood^2+c45^2-x2^2)/(2*pivotToStartWood*c45));
c45=c45*cos(theta);
c60=sqrt(pivotToStart^2+x^2-2*pivotToStart*x2*cos((180-60)*pi/180));
theta=acos((pivotToStartWood^2+c60^2-x2^2)/(2*pivotToStartWood*c60));
c60=c60*cos(theta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




FwoodTip10_0=(woodTip10_0*scaleToPivot*9.81)/(c0);
FwoodTip10_15=(woodTip10_15*scaleToPivot*9.81)/(c15);
FwoodTip10_30=(woodTip10_30*scaleToPivot*9.81)/(c30);
FwoodTip10_45=(woodTip10_45*scaleToPivot*9.81)/(c45);
FwoodTip10_60=(woodTip10_60*scaleToPivot*9.81)/(c60);

FwoodTip15_0=(woodTip15_0*scaleToPivot*9.81)/(c0);
FwoodTip15_15=(woodTip15_15*scaleToPivot*9.81)/(c15);
FwoodTip15_30=(woodTip15_30*scaleToPivot*9.81)/(c30);
FwoodTip15_45=(woodTip15_45*scaleToPivot*9.81)/(c45);
FwoodTip15_60=(woodTip15_60*scaleToPivot*9.81)/(c60);

FwoodTip20_0=(woodTip20_0*scaleToPivot*9.81)/(c0);
FwoodTip20_15=(woodTip20_15*scaleToPivot*9.81)/(c15);
FwoodTip20_30=(woodTip20_30*scaleToPivot*9.81)/(c30);
FwoodTip20_45=(woodTip20_45*scaleToPivot*9.81)/(c45);
FwoodTip20_60=(woodTip20_60*scaleToPivot*9.81)/(c60);

FwoodTip25_0=(woodTip25_0*scaleToPivot*9.81)/(c0);
FwoodTip25_15=(woodTip25_15*scaleToPivot*9.81)/(c15);
FwoodTip25_30=(woodTip25_30*scaleToPivot*9.81)/(c30);
FwoodTip25_45=(woodTip25_45*scaleToPivot*9.81)/(c45);
FwoodTip25_60=(woodTip25_60*scaleToPivot*9.81)/(c60);

FwoodTip30_0=(woodTip30_0*scaleToPivot*9.81)/(c0);
FwoodTip30_15=(woodTip30_15*scaleToPivot*9.81)/(c15);
FwoodTip30_30=(woodTip30_30*scaleToPivot*9.81)/(c30);
FwoodTip30_45=(woodTip30_45*scaleToPivot*9.81)/(c45);
FwoodTip30_60=(woodTip30_60*scaleToPivot*9.81)/(c60);