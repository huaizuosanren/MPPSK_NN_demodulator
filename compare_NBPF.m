clc;
clear all;
close all;
%Vd=[0 3 5 10];
snr=[-2 0 2 4 6 8];
%snr=[8 10 12 14 16];
Vd=0;
alpha=1.5;
M= 16 %8   %4 % 64 %  %2  %4  %2  %16   %64    %128;     %128½øÖÆ
K =24  %32  %2 % 4  % 57  % % 20  %2  %99  %198  %24    %6     %3; %2;    % every diverse is 2T
N = 400

[snrNB,peNB]=func_mppsk_demod_mppsk_chuantong (snr,M,K,N,Vd,alpha);
[snrn,pen]= func_mppsk_demod_mppsk_MP_nNBPF (snr,M,K,N,Vd,alpha);
 [snrcnn,pecnn]= func_mppsk_demod_CNN()