clear all;
close all;
clc;
snr=[-4 -2 0 2 4];
%snr=[-10 -8];
M=4
K=3
N=30
%[snr_chuantong,pe_aid3]=func_mppsk_demod_mppsk_chuantong(snr,M,K,N);
%[snr_CNN,pe_CNN]=func_mppsk_demod_CNN();
%[snr_CNN,pe_CNN3]=func_mppsk_demod_noLPF_pulse(snr,M,K,N);
n=length(snr);

% semilogy(snr_chuantong(1:n),pe_chuantong(1:n),'r--o');
%grid on;
%hold on;


M=4
K=5
N=30
[snr_CNN,pe_CNN5]=func_mppsk_demod_noLPF_pulse(snr,M,K,N);
%[snr_chuantong,pe_aid5]=func_mppsk_demod_mppsk_chuantong(snr,M,K,N);
M=4
K=7
N=30
[snr_CNN,pe_CNN7]=func_mppsk_demod_noLPF_pulse(snr,M,K,N);
%[snr_chuantong,pe_aid7]=func_mppsk_demod_mppsk_chuantong(snr,M,K,N);

%M=4
%K=15
%N=100
%[snr_chuantong,pe_aid15]=func_mppsk_demod_mppsk_chuantong(snr,M,K,N);
%M=4
%K=24
%N=100
%[snr_chuantong,pe_aid24]=func_mppsk_demod_mppsk_chuantong(snr,M,K,N);
 figure (101)
 %hold on;
 % semilogy(snr(1:n),pe_aid3(1:n),'b-*');
 % semilogy(snr(1:n),pe_CNN3(1:n),'b-*');
  hold on;
% semilogy(snr(1:n),pe_aid5(1:n),'r-o');
 semilogy(snr(1:n),pe_CNN5(1:n),'r-o');
 hold on;
% semilogy(snr(1:n),pe_aid7(1:n),'b-^');
  semilogy(snr(1:n),pe_CNN7(1:n),'b-^');
% hold on;
%  semilogy(snr(1:n),pe_aid15(1:n),'b-+');
%  hold on;
%   semilogy(snr(1:n),pe_aid24(1:n),'b-s');
 legend('K=3','K=5','K=7'); %,'K=15','K=24')
  xlabel('SNR/dB'),ylabel('SER')
%  title('diff K for AID')
   title('diff K for NN')
%  legend('AID','NN')
%legend('K=3','K=5','K=7')
grid on;
%str= sprintf('pe%d',fs/fc);
%save(str,'pe');
%str=sprintf('snr%d',fs/fc);
%save(str,'snr');
