%16 may- tip test (wood, gorilla and cheetah) temp 18.3 
%pressure 1018.29 (1017.88 to 1019.07)
clc
clear

%import the data 
run('tipWindTunnel0');
run('tipWindTunnel15');
run('tipWindTunnel30');
run('tipWindTunnel45');
run('tipWindTunnel60');
run('tipWindTunnelCalibration');

%calibrate the data
run('calibrationTipData');%remove the mean value from the rig (the top piece)

%remove the wood effect
run('calibrateTipWoodEffect2');%run script 2 and not 1 

% %calculate the force 
run('forceBlackFurTipData');
run('forceCheetahFurTipData');
run('forceWoodTip');

%calculate the drag coeff!
run('CdTipData');

%the force vectors(vector of all the forces to make plotting easier
run('dataVectors');

%run scripts to generate the graphs 
run('forcePlots');
run('CdSpeedPlots');
run('CdAnglePlots');
run('cdvsReynolds');