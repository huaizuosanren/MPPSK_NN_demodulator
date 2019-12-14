clear all;
close all;
clc;
%snr=[-6 -4 -2 0 2 4 6];
snr=[-2 0 2 4 6 8 10 12 14 16];
M=8
K=32 %3
N=400   %10 0
Vd=[0 10];
[snr_chuantong0,pe_chuantong0]=func_mppsk_demod_mppsk_chuantong(snr,M,K,N,Vd(1));
[snr_chuantong10,pe_chuantong10]=func_mppsk_demod_mppsk_chuantong(snr,M,K,N,Vd(2));
n=length(snr_chuantong0);

%[snr_CNN,pe_CNN]=func_mppsk_demod_CNN();
%[snr_CNN,pe_CNN3]=func_mppsk_demod_noLPF_pulse(snr,M,K,N);
%[snr_CNN,pe_CNN]=func_mppsk_demod_noLPF_gaussine(snr,M,K,N);
[snr_NN0,pe_NN0]= func_mppsk_demod_diffVd(Vd(1));
[snr_NN10,pe_NN10]= func_mppsk_demod_diffVd(Vd(2));

figure (102)
 semilogy(snr_chuantong0(1:n),pe_chuantong0(1:n),'r--o');
grid on;
hold on;
 semilogy(snr_chuantong10(1:n),pe_chuantong10(1:n),'r--^');
 hold on;
m=length(snr_NN0);
 semilogy(snr_NN0(1:m),pe_NN0(1:m),'b-*');
 hold on;
  semilogy(snr_NN10(1:m),pe_NN10(1:m),'b-s');

  legend('MP-Vd=0','MP-Vd=10','NN-Vd=0','NN-Vd=10')
%legend('K=3','K=5','K=7')
%grid on;
%str= sprintf('pe%d',fs/fc);
%save(str,'pe');
%str=sprintf('snr%d',fs/fc);
%save(str,'snr');
