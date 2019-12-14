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
K =[48 24 12 6] %32  %2 % 4  % 57  % % 20  %2  %99  %198  %24    %6     %3; %2;    % every diverse is 2T
N = 400
%code=[11 56 62 64]
code=[11 62 226]
M=[8 16 32 64]
[snrcnn1,pecnn1]= func_mppsk_demod_CNN(M(1),N,K(1),code(2));

[snrcnn2,pecnn2]= func_mppsk_demod_CNN(M(2),N,K(2),code(2));
[snrcnn3,pecnn3]= func_mppsk_demod_CNN(M(3),N,K(3),code(2));
[snrcnn4,pecnn4]= func_mppsk_demod_CNN(M(4),N,K(4),code(2));


%[snrcnn3,pecnn3]= func_mppsk_demod_CNN(code(4));
figure (107)
 semilogy(snrcnn1,pecnn1,'b-*');
 hold on;
 semilogy(snrcnn2,pecnn2,'b-s');
 hold on;
 semilogy(snrcnn3,pecnn3,'b-^');
 hold on;
 semilogy(snrcnn4,pecnn4,'r-^');
 legend('8','16','32','64');
  xlabel('SNR/dB'),ylabel('SER')
grid on;

[snr1,pepulse1]= func_mppsk_demod_mppsk_chuantong (snr,M(1),K(1),N,Vd,alpha);
[snr2,pepulse2]= func_mppsk_demod_mppsk_chuantong (snr,M(2),K(2),N,Vd,alpha);
[snr3,pepulse3]= func_mppsk_demod_mppsk_chuantong (snr,M(3),K(3),N,Vd,alpha);
[snr4,pepulse4]= func_mppsk_demod_mppsk_chuantong (snr,M(4),K(4),N,Vd,alpha);
figure (109)
 semilogy(snr1,pepulse1,'b-*');
 hold on;
 semilogy(snr2,pepulse2,'b-s');
 hold on;
 semilogy(snr3,pepulse3,'b-^');
 hold on;
 semilogy(snr4,pepulse4,'r-^');
 legend('8','16','32','64');
  xlabel('SNR/dB'),ylabel('SER')
grid on;