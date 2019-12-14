clear all;
close all;
clc;
%Vd=[0 3 5 10];
snr=[-6 -4 -2 0 2 4 6];
%snr=[8 10 12 14 16];
Vd=0;
alpha=1.5;
M= 16 %8   %4 % 64 %  %2  %4  %2  %16   %64    %128;     %128½øÖÆ
  %32  %2 % 4  % 57  % % 20  %2  %99  %198  %24    %6     %3; %2;    % every diverse is 2T
N = 400
%code=[11 56 62 64]
%code=[11 62 226]
code=62;
K_0 =24
[snrpnn,pepnn]= func_mppsk_demod_CNN(M,N,K_0,code);
K1=12
[snrpnSIF,pepnSIF]= func_mppsk_demod_CNN(M,N,K1,code);
K2=8
[snrMP,peMP]=func_mppsk_demod_CNN(M,N,K2,code);
K3=4
[snrel,peel]= func_mppsk_demod_CNN(M,N,K3,code);
figure (109)
 semilogy(snrpnn,pepnn,'b-*');
 hold on;
 semilogy(snrpnSIF,pepnSIF,'b-s');
 hold on;
 semilogy(snrMP,peMP,'b-^');
 hold on;
 semilogy(snrel,peel,'b-o');
 legend('K=24','K=12','K=8','K=4');
  xlabel('SNR/dB'),ylabel('SER')
grid on;
