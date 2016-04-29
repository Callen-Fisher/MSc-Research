%calculate the force on the tip of the tail covered in black fur through
%the moment arm 
%parameters
scaleToPivot=0.16;
lengthTail=0.355;
pivotToStartWood=0.12;
%fur calc
lengthFur=0.32;
a=0.009;%width of smaller edge
b=0.028;%width of bigger edge 
x2=(lengthFur/3)*(2*a+b)/(a+b);%distance to centre of area
lengthWood=0.035;
x2=x2+lengthWood;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0=sqrt(pivotToStartWood^2+x2^2-2*pivotToStartWood*x2*cos((180-0)*pi/180));%moment arm of actual force and not scale force
theta=acos((pivotToStartWood^2+c0^2-x2^2)/(2*pivotToStartWood*c0));
c0=c0*cos(theta)
c15=sqrt(pivotToStartWood^2+x2^2-2*pivotToStartWood*x2*cos((180-15)*pi/180));
theta=acos((pivotToStartWood^2+c15^2-x2^2)/(2*pivotToStartWood*c15));
c15=c15*cos(theta)
c30=sqrt(pivotToStartWood^2+x2^2-2*pivotToStartWood*x2*cos((180-30)*pi/180));
theta=acos((pivotToStartWood^2+c30^2-x2^2)/(2*pivotToStartWood*c30));
c30=c30*cos(theta)
c45=sqrt(pivotToStartWood^2+x2^2-2*pivotToStartWood*x2*cos((180-45)*pi/180));
theta=acos((pivotToStartWood^2+c45^2-x2^2)/(2*pivotToStartWood*c45));
c45=c45*cos(theta)
c60=sqrt(pivotToStartWood^2+x2^2-2*pivotToStartWood*x2*cos((180-60)*pi/180));
theta=acos((pivotToStartWood^2+c60^2-x2^2)/(2*pivotToStartWood*c60));
c60=c60*cos(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FblackFurTip10_0=(blackFurTip10_0*scaleToPivot*9.81)/(c0);
FblackFurTip10_15=(blackFurTip10_15*scaleToPivot*9.81)/(c15);
FblackFurTip10_30=(blackFurTip10_30*scaleToPivot*9.81)/(c30);
FblackFurTip10_45=(blackFurTip10_45*scaleToPivot*9.81)/(c45);
FblackFurTip10_60=(blackFurTip10_60*scaleToPivot*9.81)/(c60);

FblackFurTip15_0=(blackFurTip15_0*scaleToPivot*9.81)/(c0);
FblackFurTip15_15=(blackFurTip15_15*scaleToPivot*9.81)/(c15);
FblackFurTip15_30=(blackFurTip15_30*scaleToPivot*9.81)/(c30);
FblackFurTip15_45=(blackFurTip15_45*scaleToPivot*9.81)/(c45);
FblackFurTip15_60=(blackFurTip15_60*scaleToPivot*9.81)/(c60);

FblackFurTip20_0=(blackFurTip20_0*scaleToPivot*9.81)/(c0);
FblackFurTip20_15=(blackFurTip20_15*scaleToPivot*9.81)/(c15);
FblackFurTip20_30=(blackFurTip20_30*scaleToPivot*9.81)/(c30);
FblackFurTip20_45=(blackFurTip20_45*scaleToPivot*9.81)/(c45);
FblackFurTip20_60=(blackFurTip20_60*scaleToPivot*9.81)/(c60);

FblackFurTip25_0=(blackFurTip25_0*scaleToPivot*9.81)/(c0);
FblackFurTip25_15=(blackFurTip25_15*scaleToPivot*9.81)/(c15);
FblackFurTip25_30=(blackFurTip25_30*scaleToPivot*9.81)/(c30);
FblackFurTip25_45=(blackFurTip25_45*scaleToPivot*9.81)/(c45);
FblackFurTip25_60=(blackFurTip25_60*scaleToPivot*9.81)/(c60);

FblackFurTip30_0=(blackFurTip30_0*scaleToPivot*9.81)/(c0);
FblackFurTip30_15=(blackFurTip30_15*scaleToPivot*9.81)/(c15);
FblackFurTip30_30=(blackFurTip30_30*scaleToPivot*9.81)/(c30);
FblackFurTip30_45=(blackFurTip30_45*scaleToPivot*9.81)/(c45);
FblackFurTip30_60=(blackFurTip30_60*scaleToPivot*9.81)/(c60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FblackSprayTip10_0=(blackSprayTip10_0*scaleToPivot*9.81)/(c0);
FblackSprayTip10_15=(blackSprayTip10_15*scaleToPivot*9.81)/(c15);
FblackSprayTip10_30=(blackSprayTip10_30*scaleToPivot*9.81)/(c30);
FblackSprayTip10_45=(blackSprayTip10_45*scaleToPivot*9.81)/(c45);
FblackSprayTip10_60=(blackSprayTip10_60*scaleToPivot*9.81)/(c60);

FblackSprayTip15_0=(blackSprayTip15_0*scaleToPivot*9.81)/(c0);
FblackSprayTip15_15=(blackSprayTip15_15*scaleToPivot*9.81)/(c15);
FblackSprayTip15_30=(blackSprayTip15_30*scaleToPivot*9.81)/(c30);
FblackSprayTip15_45=(blackSprayTip15_45*scaleToPivot*9.81)/(c45);
FblackSprayTip15_60=(blackSprayTip15_60*scaleToPivot*9.81)/(c60);

FblackSprayTip20_0=(blackSprayTip20_0*scaleToPivot*9.81)/(c0);
FblackSprayTip20_15=(blackSprayTip20_15*scaleToPivot*9.81)/(c15);
FblackSprayTip20_30=(blackSprayTip20_30*scaleToPivot*9.81)/(c30);
FblackSprayTip20_45=(blackSprayTip20_45*scaleToPivot*9.81)/(c45);
FblackSprayTip20_60=(blackSprayTip20_60*scaleToPivot*9.81)/(c60);

FblackSprayTip25_0=(blackSprayTip25_0*scaleToPivot*9.81)/(c0);
FblackSprayTip25_15=(blackSprayTip25_15*scaleToPivot*9.81)/(c15);
FblackSprayTip25_30=(blackSprayTip25_30*scaleToPivot*9.81)/(c30);
FblackSprayTip25_45=(blackSprayTip25_45*scaleToPivot*9.81)/(c45);
FblackSprayTip25_60=(blackSprayTip25_60*scaleToPivot*9.81)/(c60);

FblackSprayTip30_0=(blackSprayTip30_0*scaleToPivot*9.81)/(c0);
FblackSprayTip30_15=(blackSprayTip30_15*scaleToPivot*9.81)/(c15);
FblackSprayTip30_30=(blackSprayTip30_30*scaleToPivot*9.81)/(c30);
FblackSprayTip30_45=(blackSprayTip30_45*scaleToPivot*9.81)/(c45);
FblackSprayTip30_60=(blackSprayTip30_60*scaleToPivot*9.81)/(c60);