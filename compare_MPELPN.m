clear all;
close all;
%Vd=[0 3 5 10];
snr=[-6 -4 -2 0 2 4 6];
%snr=[8 10 12 14 16];
Vd=0;
alpha=1.5;
M= 16 %8   %4 % 64 %  %2  %4  %2  %16   %64    %128;     %128½øÖÆ
K =24  %32  %2 % 4  % 57  % % 20  %2  %99  %198  %24    %6     %3; %2;    % every diverse is 2T
N = 400
%code=[11 56 62 64]
%code=[11 62 226]
code=62;
[snrpnn,pepnn]= func_mppsk_demod_CNN(M,N,K,code);
[snrpnSIF,pepnSIF]= func_mppsk_demod_CNN_nSIF(M,N,K,code);
[snrMP,peMP]=func_mppsk_demod_mppsk_chuantong (snr,M,K,N,Vd,alpha);
[snrel,peel]= func_mppsk_demod_elman(M,N,K,code);
figure (109)
 semilogy(snrpnn,pepnn,'b-*');
 hold on;
 semilogy(snrpnSIF,pepnSIF,'b-s');
 hold on;
 semilogy(snrMP,peMP,'b-^');
 hold on;
 semilogy(snrel,peel,'b-o');
 legend('SIF-FNN','FNN','SIF-MP','SIF-ELMAN');
  xlabel('SNR/dB'),ylabel('SER')
grid on;
