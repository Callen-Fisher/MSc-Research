%calculate the Cd from all the data 
a=0.009;%used to calculate the area
b=0.028;
l=0.32;
areaFur=l*(a+b)/2;

CdBlackFurTip10_0=FblackFurTip10_0(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdBlackFurTip10_15=FblackFurTip10_15(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdBlackFurTip10_30=FblackFurTip10_30(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdBlackFurTip10_45=FblackFurTip10_45(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdBlackFurTip10_60=FblackFurTip10_60(:)*2/(10*10*rho*areaFur*cos(0*pi/180));

CdBlackFurTip15_0=FblackFurTip15_0(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdBlackFurTip15_15=FblackFurTip15_15(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdBlackFurTip15_30=FblackFurTip15_30(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdBlackFurTip15_45=FblackFurTip15_45(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdBlackFurTip15_60=FblackFurTip15_60(:)*2/(15*15*rho*areaFur*cos(0*pi/180));

CdBlackFurTip20_0=FblackFurTip20_0(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdBlackFurTip20_15=FblackFurTip20_15(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdBlackFurTip20_30=FblackFurTip20_30(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdBlackFurTip20_45=FblackFurTip20_45(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdBlackFurTip20_60=FblackFurTip20_60(:)*2/(20*20*rho*areaFur*cos(0*pi/180));

CdBlackFurTip25_0=FblackFurTip25_0(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdBlackFurTip25_15=FblackFurTip25_15(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdBlackFurTip25_30=FblackFurTip25_30(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdBlackFurTip25_45=FblackFurTip25_45(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdBlackFurTip25_60=FblackFurTip25_60(:)*2/(25*25*rho*areaFur*cos(0*pi/180));

CdBlackFurTip30_0=FblackFurTip30_0(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdBlackFurTip30_15=FblackFurTip30_15(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdBlackFurTip30_30=FblackFurTip30_30(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdBlackFurTip30_45=FblackFurTip30_45(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdBlackFurTip30_60=FblackFurTip30_60(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CdBlackSprayTip10_0=FblackSprayTip10_0(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip10_15=FblackSprayTip10_15(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip10_30=FblackSprayTip10_30(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip10_45=FblackSprayTip10_45(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip10_60=FblackSprayTip10_60(:)*2/(10*10*rho*areaFur*cos(0*pi/180));

CdBlackSprayTip15_0=FblackSprayTip15_0(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip15_15=FblackSprayTip15_15(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip15_30=FblackSprayTip15_30(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip15_45=FblackSprayTip15_45(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip15_60=FblackSprayTip15_60(:)*2/(15*15*rho*areaFur*cos(0*pi/180));

CdBlackSprayTip20_0=FblackSprayTip20_0(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip20_15=FblackSprayTip20_15(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip20_30=FblackSprayTip20_30(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip20_45=FblackSprayTip20_45(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip20_60=FblackSprayTip20_60(:)*2/(20*20*rho*areaFur*cos(0*pi/180));

CdBlackSprayTip25_0=FblackSprayTip25_0(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip25_15=FblackSprayTip25_15(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip25_30=FblackSprayTip25_30(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip25_45=FblackSprayTip25_45(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip25_60=FblackSprayTip25_60(:)*2/(25*25*rho*areaFur*cos(0*pi/180));

CdBlackSprayTip30_0=FblackSprayTip30_0(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip30_15=FblackSprayTip30_15(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip30_30=FblackSprayTip30_30(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip30_45=FblackSprayTip30_45(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdBlackSprayTip30_60=FblackSprayTip30_60(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.010;
b=0.022;
l=0.24;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%should be 0.24
areaFur=l*(a+b)/2

CdCheetahSprayTip10_0=FcheetahSprayTip10_0(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip10_15=FcheetahSprayTip10_15(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip10_30=FcheetahSprayTip10_30(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip10_45=FcheetahSprayTip10_45(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip10_60=FcheetahSprayTip10_60(:)*2/(10*10*rho*areaFur*cos(0*pi/180));

CdCheetahSprayTip15_0=FcheetahSprayTip15_0(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip15_15=FcheetahSprayTip15_15(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip15_30=FcheetahSprayTip15_30(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip15_45=FcheetahSprayTip15_45(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip15_60=FcheetahSprayTip15_60(:)*2/(15*15*rho*areaFur*cos(0*pi/180));

CdCheetahSprayTip20_0=FcheetahSprayTip20_0(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip20_15=FcheetahSprayTip20_15(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip20_30=FcheetahSprayTip20_30(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip20_45=FcheetahSprayTip20_45(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip20_60=FcheetahSprayTip20_60(:)*2/(20*20*rho*areaFur*cos(0*pi/180));

CdCheetahSprayTip25_0=FcheetahSprayTip25_0(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip25_15=FcheetahSprayTip25_15(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip25_30=FcheetahSprayTip25_30(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip25_45=FcheetahSprayTip25_45(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip25_60=FcheetahSprayTip25_60(:)*2/(25*25*rho*areaFur*cos(0*pi/180));

CdCheetahSprayTip30_0=FcheetahSprayTip30_0(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip30_15=FcheetahSprayTip30_15(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip30_30=FcheetahSprayTip30_30(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip30_45=FcheetahSprayTip30_45(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdCheetahSprayTip30_60=FcheetahSprayTip30_60(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CdCheetahTip10_0=FcheetahTip10_0(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdCheetahTip10_15=FcheetahTip10_15(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdCheetahTip10_30=FcheetahTip10_30(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdCheetahTip10_45=FcheetahTip10_45(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdCheetahTip10_60=FcheetahTip10_60(:)*2/(10*10*rho*areaFur*cos(0*pi/180));

CdCheetahTip15_0=FcheetahTip15_0(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdCheetahTip15_15=FcheetahTip15_15(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdCheetahTip15_30=FcheetahTip15_30(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdCheetahTip15_45=FcheetahTip15_45(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdCheetahTip15_60=FcheetahTip15_60(:)*2/(15*15*rho*areaFur*cos(0*pi/180));

CdCheetahTip20_0=FcheetahTip20_0(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdCheetahTip20_15=FcheetahTip20_15(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdCheetahTip20_30=FcheetahTip20_30(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdCheetahTip20_45=FcheetahTip20_45(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdCheetahTip20_60=FcheetahTip20_60(:)*2/(20*20*rho*areaFur*cos(0*pi/180));

CdCheetahTip25_0=FcheetahTip25_0(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdCheetahTip25_15=FcheetahTip25_15(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdCheetahTip25_30=FcheetahTip25_30(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdCheetahTip25_45=FcheetahTip25_45(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdCheetahTip25_60=FcheetahTip25_60(:)*2/(25*25*rho*areaFur*cos(0*pi/180));

CdCheetahTip30_0=FcheetahTip30_0(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdCheetahTip30_15=FcheetahTip30_15(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdCheetahTip30_30=FcheetahTip30_30(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdCheetahTip30_45=FcheetahTip30_45(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdCheetahTip30_60=FcheetahTip30_60(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.006;
b=0.024;
l=0.355;
areaFur=l*(a+b)/2;

CdWoodTip10_0=FwoodTip10_0(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdWoodTip10_15=FwoodTip10_15(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdWoodTip10_30=FwoodTip10_30(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdWoodTip10_45=FwoodTip10_45(:)*2/(10*10*rho*areaFur*cos(0*pi/180));
CdWoodTip10_60=FwoodTip10_60(:)*2/(10*10*rho*areaFur*cos(0*pi/180));

CdWoodTip15_0=FwoodTip15_0(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdWoodTip15_15=FwoodTip15_15(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdWoodTip15_30=FwoodTip15_30(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdWoodTip15_45=FwoodTip15_45(:)*2/(15*15*rho*areaFur*cos(0*pi/180));
CdWoodTip15_60=FwoodTip15_60(:)*2/(15*15*rho*areaFur*cos(0*pi/180));

CdWoodTip20_0=FwoodTip20_0(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdWoodTip20_15=FwoodTip20_15(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdWoodTip20_30=FwoodTip20_30(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdWoodTip20_45=FwoodTip20_45(:)*2/(20*20*rho*areaFur*cos(0*pi/180));
CdWoodTip20_60=FwoodTip20_60(:)*2/(20*20*rho*areaFur*cos(0*pi/180));

CdWoodTip25_0=FwoodTip25_0(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdWoodTip25_15=FwoodTip25_15(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdWoodTip25_30=FwoodTip25_30(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdWoodTip25_45=FwoodTip25_45(:)*2/(25*25*rho*areaFur*cos(0*pi/180));
CdWoodTip25_60=FwoodTip25_60(:)*2/(25*25*rho*areaFur*cos(0*pi/180));

CdWoodTip30_0=FwoodTip30_0(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdWoodTip30_15=FwoodTip30_15(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdWoodTip30_30=FwoodTip30_30(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdWoodTip30_45=FwoodTip30_45(:)*2/(30*30*rho*areaFur*cos(0*pi/180));
CdWoodTip30_60=FwoodTip30_60(:)*2/(30*30*rho*areaFur*cos(0*pi/180));