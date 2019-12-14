clc;
clear all;
close all;
%***************************************MPPSK***********************
M=2 %8   %64  %2 %2   %3 %2  %64    %128;     %128jizhi
K =200 %  1; % %=57 % 48  %6  %8  %8  %133  %198  %6     %3; %2;    % every diverse is 2T
N =400  %399 %400  %20   %400;   %260;   % every base band signal is 260T
A=1;     % amplifier is 1
B=1;
rg=0;   
thita=pi;  % the diverse angle is pi
fc=20e3;     
fs=fc*8;   
n1=N*fs/fc;    % fs/fc is the number of samples for one timeperiod of carrier. N*fs/fc is the number of samples for one time period for base band 
n2=(1-rg)*K*fs/fc; 

   t0 = (0:1:n1 - 1)/fs;
    f0 = A*sin(2*pi*fc*t0);    
    t1 = (0:1:n2 - 1)/fs;
    f1 = B*sin(2*pi*fc*t1 + thita);

Nsymbols =20;
    symbol_sample_stream_num=Nsymbols*n1;
        mppsk=zeros(Nsymbols,n1);
        s_tx = zeros(1,Nsymbols);
        s_tx = randi([0 M - 1],1,Nsymbols);
    
 %% ******** get the modemed signals****************************
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
 
 %% *********************SIF*****************       
           b =[1  -1.414423081535762   1];     %fs/fc=8
     a =[1  -3.987007081485477   7.949450640987077  -9.385437352983406   7.009070118114206  -3.099682963193266   0.685515443313961];
          zf2 = zeros(max(length(b),length(a))-1,1);
            [mppskfilter,zf2]=filter(b,a,mppsk,zf2);  
        figure(100);
        subplot(2,1,1);
        plot(mppsk);
        subplot(2,1,2);
        plot(mppskfilter);