function [ F ] = partialNumericDiff(states )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



F= [[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    f(states,0.01);
    [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
    [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]];





% d_tp1=0.01;
% d_tp2=0.01;
% d_ty1=0.01;
% d_ty2=0.01;
% d_dsi1=0.0001;
% d_dth1=0.01;
% d_dph1=0.0001;
% d_dsi2=0.0001;
% d_dth2=0.01;
% d_dph2=0.0001;
% d_dth3=0.01;
% d_dph3=0.0001;
% d_dth4=0.01;
% d_dph4=0.0001;
% d_si1=0.0001;
% d_th1=0.01;
% d_ph1=0.0001;
% d_si2=0.0001;
% d_th2=0.01;
% d_ph2=0.0001;
% d_th3=0.01;
% d_ph3=0.0001;
% d_th4=0.01;
% d_ph4=0.0001;
% 
% F=[(current-previous)/d_tp1;
%    (current-previous)/d_ty1;
%    (current-previous)/d_tp2;
%    (current-previous)/d_ty2;
%    (current-previous)/d_dsi1;
%    (current-previous)/d_dth1;
%    (current-previous)/d_dph1;
%    (current-previous)/d_dsi2;
%    (current-previous)/d_dth2;
%    (current-previous)/d_dph2;
%    (current-previous)/d_dth3;
%    (current-previous)/d_dph3;
%    (current-previous)/d_dth4;
%    (current-previous)/d_dph4;
%    (current-previous)/d_si1;
%    (current-previous)/d_th1;
%    (current-previous)/d_ph1;
%    (current-previous)/d_si2;
%    (current-previous)/d_th2;
%    (current-previous)/d_ph2;
%    (current-previous)/d_th3;
%    (current-previous)/d_ph3;
%    (current-previous)/d_th4;
%    (current-previous)/d_ph4];

end

