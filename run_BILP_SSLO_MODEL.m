%clear
clc
%load('BILP_Grade.mat')
%load('BILP_EBV.mat')
% 
% EBV = -ones(10,10);
% EBV(3:10,3:8) = 1;
% EBV(7:10,7:10) = 1;
% % EBV(12:17,3:10) = 1;
% % EBV(4:11,9:16) = 1;
% Grade = zeros(10,10);
% Grade(3:10,3:8) = 0.5 + rand(8,6)/2;
% Grade(7:10,7:10) = 0.5 + rand(4,4)/2;
% % Grade(12:17,3:10) = 0.6 + rand(6,8)/2;
% % Grade(4:11,9:16) = 0.6 + rand(8,8)/2;
Cutoff = 3;
[sol,xFlag,OFV,Tm,EBV,Grade,stp,mm,mp] = BILP_SSLO(EBV,Grade,Cutoff);