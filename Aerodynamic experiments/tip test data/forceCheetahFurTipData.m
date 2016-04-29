%calculate the force on the tail covered in cheetah fur through the
%moment arm 
%parameters
scaleToPivot=0.16;
pivotToStartWood=0.12;
lengthTail=0.355;
%fur calc
lengthFur=0.30;
a=0.010;%width of smaller edge
b=0.022;%width of bigger edge 
x2=(lengthFur/3)*(2*a+b)/(a+b);%distance to centre of area
lengthWood=0.1;
x2=x2+lengthWood;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0=sqrt(pivotToStartWood^2+x2^2-2*pivotToStartWood*x2*cos((180-0)*pi/180));%moment arm of actual force and not scale force
theta=acos((pivotToStartWood^2+c0^2-x2^2)/(2*pivotToStartWood*c0));%angle to make fdrag perpendicular
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




FcheetahSprayTip10_0=(cheetahSprayTip10_0*scaleToPivot*9.81)/(c0);
FcheetahSprayTip10_15=(cheetahSprayTip10_15*scaleToPivot*9.81)/(c15);
FcheetahSprayTip10_30=(cheetahSprayTip10_30*scaleToPivot*9.81)/(c30);
FcheetahSprayTip10_45=(cheetahSprayTip10_45*scaleToPivot*9.81)/(c45);
FcheetahSprayTip10_60=(cheetahSprayTip10_60*scaleToPivot*9.81)/(c60);

FcheetahSprayTip15_0=(cheetahSprayTip15_0*scaleToPivot*9.81)/(c0);
FcheetahSprayTip15_15=(cheetahSprayTip15_15*scaleToPivot*9.81)/(c15);
FcheetahSprayTip15_30=(cheetahSprayTip15_30*scaleToPivot*9.81)/(c30);
FcheetahSprayTip15_45=(cheetahSprayTip15_45*scaleToPivot*9.81)/(c45);
FcheetahSprayTip15_60=(cheetahSprayTip15_60*scaleToPivot*9.81)/(c60);

FcheetahSprayTip20_0=(cheetahSprayTip20_0*scaleToPivot*9.81)/(c0);
FcheetahSprayTip20_15=(cheetahSprayTip20_15*scaleToPivot*9.81)/(c15);
FcheetahSprayTip20_30=(cheetahSprayTip20_30*scaleToPivot*9.81)/(c30);
FcheetahSprayTip20_45=(cheetahSprayTip20_45*scaleToPivot*9.81)/(c45);
FcheetahSprayTip20_60=(cheetahSprayTip20_60*scaleToPivot*9.81)/(c60);

FcheetahSprayTip25_0=(cheetahSprayTip25_0*scaleToPivot*9.81)/(c0);
FcheetahSprayTip25_15=(cheetahSprayTip25_15*scaleToPivot*9.81)/(c15);
FcheetahSprayTip25_30=(cheetahSprayTip25_30*scaleToPivot*9.81)/(c30);
FcheetahSprayTip25_45=(cheetahSprayTip25_45*scaleToPivot*9.81)/(c45);
FcheetahSprayTip25_60=(cheetahSprayTip25_60*scaleToPivot*9.81)/(c60);

FcheetahSprayTip30_0=(cheetahSprayTip30_0*scaleToPivot*9.81)/(c0);
FcheetahSprayTip30_15=(cheetahSprayTip30_15*scaleToPivot*9.81)/(c15);
FcheetahSprayTip30_30=(cheetahSprayTip30_30*scaleToPivot*9.81)/(c30);
FcheetahSprayTip30_45=(cheetahSprayTip30_45*scaleToPivot*9.81)/(c45);
FcheetahSprayTip30_60=(cheetahSprayTip30_60*scaleToPivot*9.81)/(c60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FcheetahTip10_0=(cheetahTip10_0*scaleToPivot*9.81)/(c0);
FcheetahTip10_15=(cheetahTip10_15*scaleToPivot*9.81)/(c15);
FcheetahTip10_30=(cheetahTip10_30*scaleToPivot*9.81)/(c30);
FcheetahTip10_45=(cheetahTip10_45*scaleToPivot*9.81)/(c45);
FcheetahTip10_60=(cheetahTip10_60*scaleToPivot*9.81)/(c60);

FcheetahTip15_0=(cheetahTip15_0*scaleToPivot*9.81)/(c0);
FcheetahTip15_15=(cheetahTip15_15*scaleToPivot*9.81)/(c15);
FcheetahTip15_30=(cheetahTip15_30*scaleToPivot*9.81)/(c30);
FcheetahTip15_45=(cheetahTip15_45*scaleToPivot*9.81)/(c45);
FcheetahTip15_60=(cheetahTip15_60*scaleToPivot*9.81)/(c60);

FcheetahTip20_0=(cheetahTip20_0*scaleToPivot*9.81)/(c0);
FcheetahTip20_15=(cheetahTip20_15*scaleToPivot*9.81)/(c15);
FcheetahTip20_30=(cheetahTip20_30*scaleToPivot*9.81)/(c30);
FcheetahTip20_45=(cheetahTip20_45*scaleToPivot*9.81)/(c45);
FcheetahTip20_60=(cheetahTip20_60*scaleToPivot*9.81)/(c60);

FcheetahTip25_0=(cheetahTip25_0*scaleToPivot*9.81)/(c0);
FcheetahTip25_15=(cheetahTip25_15*scaleToPivot*9.81)/(c15);
FcheetahTip25_30=(cheetahTip25_30*scaleToPivot*9.81)/(c30);
FcheetahTip25_45=(cheetahTip25_45*scaleToPivot*9.81)/(c45);
FcheetahTip25_60=(cheetahTip25_60*scaleToPivot*9.81)/(c60);

FcheetahTip30_0=(cheetahTip30_0*scaleToPivot*9.81)/(c0);
FcheetahTip30_15=(cheetahTip30_15*scaleToPivot*9.81)/(c15);
FcheetahTip30_30=(cheetahTip30_30*scaleToPivot*9.81)/(c30);
FcheetahTip30_45=(cheetahTip30_45*scaleToPivot*9.81)/(c45);
FcheetahTip30_60=(cheetahTip30_60*scaleToPivot*9.81)/(c60);