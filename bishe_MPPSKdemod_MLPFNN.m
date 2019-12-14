%clc;
%clear all;
%close all;
%Vd=[0 3 5 10];
%snr=[-2 0 2 4 6 8];
%snr=[8 10 12 14 16];
snr=[-6 -4 -2 0 2 4 6];
Vd=0;
alpha=1.5;
%M= 16 %8   %4 % 64 %  %2  %4  %2  %16   %64    %128;     %128½øÖÆ
N = 400
M=16;
K=24;
Vd = 0;
code = 62;
[snrcnn1,pecnn1]= func_mppsk_demod_CNN(M,N,K,code);
[snr_chuantong10,pe_chuantong10]=func_mppsk_demod_mppsk_chuantong(snr,M,K,N,Vd);
[snr,pecnn2]= func_mppsk_demod_CNN_nSIF(M,N,K,code)
figure (106)
 semilogy(snrcnn1,pecnn1,'b-*');
 hold on;
 semilogy(snr_chuantong10,pe_chuantong10,'b-^');
 hold on;
 semilogy(snr,pecnn2,'b-s');
 legend('proposed demodulation scheme','demodulation scheme in [7]','proposed demodulation scheme without SIF');
  xlabel('SNR/dB'),ylabel('SER')
grid on;

figure (107)
 semilogy(snrcnn1,pecnn1,'b-*');
 legend('16PPSK');
 xlabel('SNR/dB'),ylabel('SER')
grid on;
