clc;
clear all;
close all;
%Vd=[0 3 5 10]
Vd=0;
alpha=[1.2 1.8];
snr=[-2 0 2 4 6 8];
M= 16 %8   %4 % 64 %  %2  %4  %2  %16   %64    %128;     %128½øÖÆ
K = 24 %32  %2 % 4  % 57  % % 20  %2  %99  %198  %24    %6     %3; %2;    % every diverse is 2T
N = 400 %10 %256   %399 %  %180  %10  %400;   %260;   % every base band signal is 260T
code=62;
[snr0,pe0]= func_mppsk_demod_diffVd(Vd,alpha(1),code);
[snr3,pe3]= func_mppsk_demod_diffVd(Vd,alpha(2),code);
%[snr5,pe5]= func_mppsk_demod_diffVd(Vd,alpha(3),10);
%[snr10,pe10]= func_mppsk_demod_diffVd(Vd,alpha(4),10);
[snr5,pe5]=func_mppsk_demod_mppsk_chuantong (snr,M,K,N,Vd,alpha(1));
[snr10,pe10]= func_mppsk_demod_mppsk_chuantong (snr,M,K,N,Vd,alpha(2));

%[snr0,pe0]= func_mppsk_demod_diffVd(Vd,alpha(1),10);
%[snr3,pe3]= func_mppsk_demod_diffVd(Vd,alpha(1),11);
%[snr5,pe5]= func_mppsk_demod_diffVd(Vd,alpha(1),64);
%[snr10,pe10]= func_mppsk_demod_diffVd(Vd,alpha(1),128);
figure(105)
 semilogy(snr0,pe0,'b-*');
 hold on;
 semilogy(snr3,pe3,'b-o');
 hold on;
 semilogy(snr5,pe5,'b-s');
 hold on;
 semilogy(snr10,pe10,'b->');
%  legend('Vd=0','Vd=3','Vd=5','Vd=10');
   legend('FNN-alpha=1.2','alpha=1.6','alpha=1.7','alpha=1.8');
 legend('FNN-alpha=1.2','FNN-alpha=1.8','MP-alpha=1.2','MP-alpha=1.8');
  xlabel('SNR/dB'),ylabel('SER')
grid on;