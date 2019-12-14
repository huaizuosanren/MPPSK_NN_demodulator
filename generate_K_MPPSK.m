clc;
clear all;
close all;
%% **********generate MPPSK********************
M=16 %8   %64  %2 %2   %3 %2  %64    %128;     %128jizhi
K =24  %=57 % 48  %6  %8  %8  %133  %198  %6     %3; %2;    % every diverse is 2T
N =400  %399 %400  %20   %400;   %260;   % every base band signal is 260T
A=1;     % amplifier is 1
B=1;
rg=0;   
thita=pi;  % the diverse angle is pi
fc=20e3;     
fs=fc*8;   
n1=N*fs/fc;    % fs/fc is the number of samples for one timeperiod of carrier. N*fs/fc is the number of samples for one time period for base band 
n2=(1-rg)*K*fs/fc;   %the number of samples for one diverse time period.
mm=[];
ff=[];
p_signal=0.5;
%snr=[-12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5];
%snr=[0 1 2 3 4 5 6 7 8];
%snr=[-5 -10 -8 -6 -4 -2 ];
%snr=[-5 -10 -8 -6 -4 -2 0];
%snr=[-6 -5 -4 -3 -2 -1 0 1 2];
%Nblocks = [40 40 80 80 80 160 160 160 160 200 20 40 40 40 80 80 80 160 160 160 160 160 160];
Nblocks = [20 20 20 20 20 20 20 20 20 20 20 20 20 20 20];
%Nblocks = [20 100 500 1000 1000 1000 2000 5000 5000 5000];

pe=[];
mppskfilter=[];

Nsymbols =5000;

 % b=[1,-1.9021496572560159 ,1];     % fs/fc=20
 %   a=[1,-3.6512241163814698 , 5.1727286626648894 ,-3.3577825365961242, 0.84572301542400019];
 
   b =[1  -1.414423081535762   1];     %fs/fc=8
     a =[1  -3.987007081485477   7.949450640987077  -9.385437352983406   7.009070118114206  -3.099682963193266   0.685515443313961];
    % Hd=lowfilter_eb_8;
   %  b_env=conv( conv( Hd.sosMatrix(1,1:3),Hd.sosMatrix(2,1:3) ),Hd.sosMatrix(3,1:3) )*Hd.ScaleValues(1)*Hd.ScaleValues(2)*Hd.ScaleValues(3)*Hd.ScaleValues(4);
   %  a_env=conv( conv( Hd.sosMatrix(1,4:6),Hd.sosMatrix(2,4:6) ) ,Hd.sosMatrix(3,4:6));
 
 
  %    load SOS;
  %  load sosMatrix;
    load num_low;
    Hd=tianxian;
%    [b_shape,a_shape]=sos2tf(SOS,G);
 %  [b_shape,a_shape]=sos2tf(sosMatrix)  %,G);
    fprintf('have not been filtered  sample ratio£º8\n \n');
    
    
%% initialize the variables---------------------------------------------------------------------
        symbol_sample_stream_num=Nsymbols*n1;
        mppsk=zeros(Nsymbols,n1);
        s_tx = zeros(1,Nsymbols);  %,[0 M - 1]);

        
    
    t0 = (0:1:n1 - 1)/fs;
    f0 = A*sin(2*pi*fc*t0);    
    t1 = (0:1:n2 - 1)/fs;
    f1 = B*sin(2*pi*fc*t1 + thita);
    
  %  err = 0;
    s_rcv=[];
    s_tx_total=[];
    zf1 = zeros(max(length(b),length(a))-1,1);
    zf2 = zeros(max(length(b),length(a))-1,1);
    zf_env1 = zeros(max(length(Num_low),1)-1,1);
    zf_env2 = zeros(max(length(Num_low),1)-1,1);
  %  zf_shape=zeros(max(length(b_shape),length(a_shape))-1,1);
 
        symbol_sample_stream_num=Nsymbols*n1;
        mppsk=zeros(Nsymbols,n1);
        %s_tx = randint(1,Nsymbols,[0 M - 1]);
        s_tx = randi([0 M - 1],1,Nsymbols);
    
 %% ******** get the modemed signals and add noise to it****************************
        for num = 1:Nsymbols
            if( s_tx(num) == 0)
                mppsk(num, :) = f0;                      
            else                                        
                k = s_tx(num);                    
                head = f0(1:(k - 1)*K*(fs/fc));        
                middle=f1;
                tail = f0(1:(N - (k - rg)*K)*(fs/fc));   
                mppsk(num, :) = [head middle tail];     
            end
        end
        mppsk = reshape(mppsk', 1, symbol_sample_stream_num);
  %*************pass through the tianxian*************************      
        %[mppsk,zf_shape]=filter(b_shape,a_shape,mppsk);
        NB_mppsk=filter(Hd,mppsk);
%**************pass the SIF*************        
       [mppskfilter,zf2]=filter(b,a,mppsk,zf2);      
       [NB_mppskfilter,zf2]=filter(b,a,NB_mppsk,zf2);   


 %% *****************plot the fiure*************************** 
  %--- Initialization ---------------------------------------------------
    figure(100);
    
    
    clf(100);
    
    timeScale = 0 : 1/fs : 0.1;    
    
    %--- Time domain plot -------------------------------------------------
    subplot(2, 2, 1);
   % plot(1000 * timeScale(1:round(samplesPerCode/50)), ...
    %     data(1:round(samplesPerCode/50)));
      plot(1000* timeScale, ...
         mppsk(1:length(timeScale)));
     
     
    axis tight;
    grid on;
    title ('MPPSK Time domain plot');
    xlabel('Time (ms)'); ylabel('Amplitude');
    %-----------after SIF time domin---------
        subplot(2, 2, 2);
   % plot(1000 * timeScale(1:round(samplesPerCode/50)), ...
    %     data(1:round(samplesPerCode/50)));
      plot(1000* timeScale, ...
         mppskfilter(1:length(timeScale)));
     
     
    axis tight;
    grid on;
    title ('MPPSK after SIF');
    xlabel('Time (ms)'); ylabel('Amplitude');
    
   %---NBPF Time domain plot -------------------------------------------------
    subplot(2, 2, 3);
   % plot(1000 * timeScale(1:round(samplesPerCode/50)), ...
    %     data(1:round(samplesPerCode/50)));
      plot(1000* timeScale, ...
         NB_mppsk(1:length(timeScale)));
     
     
    axis tight;
    grid on;
    title ('NBL MPPSK Time domain plot');
    xlabel('Time (ms)'); ylabel('Amplitude');
    %-----------NBPF after SIF time domin---------
        subplot(2, 2, 4);
   % plot(1000 * timeScale(1:round(samplesPerCode/50)), ...
    %     data(1:round(samplesPerCode/50)));
      plot(1000* timeScale, ...
         NB_mppskfilter(1:length(timeScale)));
     
     
    axis tight;
    grid on;
    title ('NBL MPPSK after SIF');
    xlabel('Time (ms)'); ylabel('Amplitude');
       
    
    %--- Frequency domain plot --------------------------------------------
    
      figure(66)
   
   % [pwm,fm]=pwelch(mppsk(1:16)-mean(mppsk(1:16)),fs);   %, 16384, 1024, 2048, settings.samplingFreq/1e6)
  % pwm=pwelch(mppsk(1:16)-mean(mppsk(1:16)));
  % plot(10*log10(pwm),'b');
   pwelch(mppsk-mean(mppsk));
    axis tight;
    grid on;
    title ('Frequency domain plot');
    xlabel('normalized Frequency'); ylabel('Magnitude');
    hold on;
   % [pwn,fn]=pwelch(NB_mppsk(1:16)-mean(NB_mppsk(1:16)),fs);
   % plot(fn,pwn,'r');
     %pwn=pwelch(NB_mppsk(1:16)-mean(NB_mppsk(1:16)));
     pwelch(NB_mppsk-mean(NB_mppsk));
  % plot(10*log10(pwn),'r');
    %legend('MPPSK','NBL MPPSK');
  
    
    figure(101)
   
   % [pwm,fm]=pwelch(mppsk(1:16)-mean(mppsk(1:16)),fs);   %, 16384, 1024, 2048, settings.samplingFreq/1e6)
  % pwm=pwelch(mppsk(1:16)-mean(mppsk(1:16)));
  % plot(10*log10(pwm),'b');
   pwelch(mppsk-mean(mppsk));
    axis tight;
    grid on;
    title ('Frequency domain plot');
    xlabel('normalized Frequency'); ylabel('Magnitude');
    hold on;
   % [pwn,fn]=pwelch(NB_mppsk(1:16)-mean(NB_mppsk(1:16)),fs);
   % plot(fn,pwn,'r');
     %pwn=pwelch(NB_mppsk(1:16)-mean(NB_mppsk(1:16)));
     pwelch(NB_mppsk-mean(NB_mppsk));
  % plot(10*log10(pwn),'r');
    %legend('MPPSK','NBL MPPSK');
    hold on;
 
    %% ************K=8****************************
    clear all;
    M=16 %8   %64  %2 %2   %3 %2  %64    %128;     %128jizhi
K =8  %24  %=57 % 48  %6  %8  %8  %133  %198  %6     %3; %2;    % every diverse is 2T
N =400  %399 %400  %20   %400;   %260;   % every base band signal is 260T
A=1;     % amplifier is 1
B=1;
rg=0;   
thita=pi;  % the diverse angle is pi
fc=20e3;     
fs=fc*8;   
n1=N*fs/fc;    % fs/fc is the number of samples for one timeperiod of carrier. N*fs/fc is the number of samples for one time period for base band 
n2=(1-rg)*K*fs/fc;   %the number of samples for one diverse time period.
mm=[];
ff=[];
p_signal=0.5;
%snr=[-12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5];
%snr=[0 1 2 3 4 5 6 7 8];
%snr=[-5 -10 -8 -6 -4 -2 ];
%snr=[-5 -10 -8 -6 -4 -2 0];
%snr=[-6 -5 -4 -3 -2 -1 0 1 2];
%Nblocks = [40 40 80 80 80 160 160 160 160 200 20 40 40 40 80 80 80 160 160 160 160 160 160];
Nblocks = [20 20 20 20 20 20 20 20 20 20 20 20 20 20 20];
%Nblocks = [20 100 500 1000 1000 1000 2000 5000 5000 5000];

pe=[];
mppskfilter=[];

Nsymbols =5000;

 % b=[1,-1.9021496572560159 ,1];     % fs/fc=20
 %   a=[1,-3.6512241163814698 , 5.1727286626648894 ,-3.3577825365961242, 0.84572301542400019];
 
   b =[1  -1.414423081535762   1];     %fs/fc=8
     a =[1  -3.987007081485477   7.949450640987077  -9.385437352983406   7.009070118114206  -3.099682963193266   0.685515443313961];
    % Hd=lowfilter_eb_8;
   %  b_env=conv( conv( Hd.sosMatrix(1,1:3),Hd.sosMatrix(2,1:3) ),Hd.sosMatrix(3,1:3) )*Hd.ScaleValues(1)*Hd.ScaleValues(2)*Hd.ScaleValues(3)*Hd.ScaleValues(4);
   %  a_env=conv( conv( Hd.sosMatrix(1,4:6),Hd.sosMatrix(2,4:6) ) ,Hd.sosMatrix(3,4:6));
 
 
  %    load SOS;
  %  load sosMatrix;
    load num_low;
    Hd=tianxian;
%    [b_shape,a_shape]=sos2tf(SOS,G);
 %  [b_shape,a_shape]=sos2tf(sosMatrix)  %,G);
    fprintf('have not been filtered  sample ratio£º8\n \n');
    
    
%% initialize the variables---------------------------------------------------------------------
        symbol_sample_stream_num=Nsymbols*n1;
        mppsk=zeros(Nsymbols,n1);
        s_tx = zeros(1,Nsymbols);  %,[0 M - 1]);

        
    
    t0 = (0:1:n1 - 1)/fs;
    f0 = A*sin(2*pi*fc*t0);    
    t1 = (0:1:n2 - 1)/fs;
    f1 = B*sin(2*pi*fc*t1 + thita);
    
  %  err = 0;
    s_rcv=[];
    s_tx_total=[];
    zf1 = zeros(max(length(b),length(a))-1,1);
    zf2 = zeros(max(length(b),length(a))-1,1);
    zf_env1 = zeros(max(length(Num_low),1)-1,1);
    zf_env2 = zeros(max(length(Num_low),1)-1,1);
  %  zf_shape=zeros(max(length(b_shape),length(a_shape))-1,1);
 
        symbol_sample_stream_num=Nsymbols*n1;
        mppsk=zeros(Nsymbols,n1);
        %s_tx = randint(1,Nsymbols,[0 M - 1]);
        s_tx = randi([0 M - 1],1,Nsymbols);
    
 %% ******** get the modemed signals and add noise to it****************************
        for num = 1:Nsymbols
            if( s_tx(num) == 0)
                mppsk(num, :) = f0;                      
            else                                        
                k = s_tx(num);                    
                head = f0(1:(k - 1)*K*(fs/fc));        
                middle=f1;
                tail = f0(1:(N - (k - rg)*K)*(fs/fc));   
                mppsk(num, :) = [head middle tail];     
            end
        end
        mppsk = reshape(mppsk', 1, symbol_sample_stream_num);
  %*************pass through the tianxian*************************      
        %[mppsk,zf_shape]=filter(b_shape,a_shape,mppsk);
        NB_mppsk=filter(Hd,mppsk);
%**************pass the SIF*************        
       [mppskfilter,zf2]=filter(b,a,mppsk,zf2);      
       [NB_mppskfilter,zf2]=filter(b,a,NB_mppsk,zf2);   
%% **************plo***************************
   pwelch(mppsk-mean(mppsk));
    axis tight;
    grid on;
    title ('Frequency domain plot');
    xlabel('normalized Frequency'); ylabel('Magnitude');
    hold on;
   % [pwn,fn]=pwelch(NB_mppsk(1:16)-mean(NB_mppsk(1:16)),fs);
   % plot(fn,pwn,'r');
     %pwn=pwelch(NB_mppsk(1:16)-mean(NB_mppsk(1:16)));
     pwelch(NB_mppsk-mean(NB_mppsk));
  % plot(10*log10(pwn),'r');
    %legend('MPPSK','NBL MPPSK');
    hold on;
    
    figure(67)
       pwelch(mppsk-mean(mppsk));
    axis tight;
    grid on;
    title ('Frequency domain plot');
    xlabel('normalized Frequency'); ylabel('Magnitude');
    hold on;
   % [pwn,fn]=pwelch(NB_mppsk(1:16)-mean(NB_mppsk(1:16)),fs);
   % plot(fn,pwn,'r');
     %pwn=pwelch(NB_mppsk(1:16)-mean(NB_mppsk(1:16)));
     pwelch(NB_mppsk-mean(NB_mppsk));