clc;
clear all;
close all;
%Vd=[0 3 5 10];
%snr=[-2 0 2 4 6 8];
snr=[8 10 12 14 16];
Vd=0;
alpha=[1.5 1.6 1.7 1.8]
M=8   %4 % 64 %  %2  %4  %2  %16   %64    %128;     %128½øÖÆ
K =32  %2 % 4  % 57  % % 20  %2  %99  %198  %24    %6     %3; %2;    % every diverse is 2T
N = 400

[snr0,pe0]=func_mppsk_demod_mppsk_chuantong (snr,M,K,N,Vd,alpha(1));
[snr3,pe3]= func_mppsk_demod_mppsk_chuantong (snr,M,K,N,Vd,alpha(2));
[snr5,pe5]=func_mppsk_demod_mppsk_chuantong (snr,M,K,N,Vd,alpha(3));
[snr10,pe10]=func_mppsk_demod_mppsk_chuantong (snr,M,K,N,Vd,alpha(4));
figure(105)
 semilogy(snr0,pe0,'b-*');
 hold on;
 semilogy(snr3,pe3,'b-o');
 hold on;
 semilogy(snr5,pe5,'b-s');
 hold on;
 semilogy(snr10,pe10,'b->');
%  legend('Vd=0','Vd=3','Vd=5','Vd=10');
  legend('alpha=1.5','alpha=1.6','alpha=1.7','alpha=1.8');

  xlabel('SNR/dB'),ylabel('SER')
 % title('different Vd for a-stable pulse noise');
grid on;