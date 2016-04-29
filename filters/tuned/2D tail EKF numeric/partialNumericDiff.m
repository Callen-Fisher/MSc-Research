function [ F ] = partialNumericDiff(states )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

F= [f1(states,0.01);
    [1 0 0 0];
    [0 1 0 0]];

% 
% d_dth=0.01;
% d_th=0.01;
% 
% F=[(current-previous)/d_dth;(current-previous)/d_dth;(current-previous)/d_th;(current-previous)/d_th];

end

