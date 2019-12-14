clc;
clear all;
close all;
%Vd=[0 3 5 10];
%snr=[-2 0 2 4 6 8];
%snr=[8 10 12 14 16];
snr=[-6 -4 -2 0 2 4 6];
Vd=0;
alpha=1.5;
%M= 16 %8   %4 % 64 %  %2  %4  %2  %16   %64    %128;     %128½øÖÆ
K =[24 24 24 6] %32  %2 % 4  % 57  % % 20  %2  %99  %198  %24    %6     %3; %2;    % every diverse is 2T
N = 400
%code=[11 56 62 64]
code=[11 62 226]
M=[8 4 32 64]
[snrcnn1,pecnn1]= func_mppsk_demod_CNN(M(1),N,K(1),code(2));
[snrcnn2,pecnn2]= func_mppsk_demod_CNN(M(2),N,K(2),code(2));
figure (107)
 semilogy(snrcnn1,pecnn1,'b-*');
 semilogy(snrcnn2,pecnn2,'b-^');
  xlabel('SNR/dB'),ylabel('SER')
grid on;