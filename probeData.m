function [data,mppsk,s_tx,k]= probeData(varargin)
%Function plots raw data information: time domain plot, a frequency domain
%plot and a histogram. 
%
%The function can be called in two ways:
%   probeData(settings)
% or
%   probeData(fileName, settings)
%
%   Inputs:
%       fileName        - name of the data file. File name is read from
%                       settings if parameter fileName is not provided.
%
%       settings        - receiver settings. Type of data file, sampling
%                       frequency and the default filename are specified
%                       here. 

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

% CVS record:
% $Id: probeData.m,v 1.1.2.7 2006/08/22 13:46:00 dpl Exp $

    
%% Generate plot of raw data ==============================================

if (nargin == 1)
    settings = deal(varargin{1});

else
    error('Incorect number of arguments');
end

    samplesPerCode = round(settings.samplingFreq / ...
                           (1e3));
      
                        % Read data for acquisition. 11ms of signal are needed for the
        % fine   麓脫虏露禄帽脰脨露脕脠隆脢媒戮脻11ms
        % frequency estimation
        %麓脫脦脛录镁脰脨露脕脠隆脨脜潞脜 data = fread(fid, 11*samplesPerCode, settings.dataType)';
        %Generate the carrier frequency 脫脙碌楼脭脴虏篓脨脜潞脜脠隆麓煤脭颅脦脛录镁脨脜潞脜 -----------
           
      %%Prepare for generating signals  
     
            time=settings.msToProcess;  %1000ms   
    carrFreq_rcv= 20e3; %15.580030e6;
            blksize_rcv=time*samplesPerCode;  %about 1000ms's signal,there are 1000 c/a periods 50 base periods
            time_rcv    = (0:blksize_rcv) ./ settings.samplingFreq;
           % caCodeForSamples=zeros(1,samplesPerCode);
 %% initial the preparation for MPPSK     
 M=64    %128;     %128进制
K =6     %3; %2;    % every diverse is 2T
N =400;   %260;   % every base band signal is 260T
A=1;     % amplifier is 1
B=1;
rg=0;   
thita=pi;  % the diverse angle is pi


fc=carrFreq_rcv;     
fs=fc*8;   
n1=N*fs/fc;    % fs/fc is the number of samples for one timeperiod of carrier. N*fs/fc is the number of samples for one time period for base band 
n2=(1-rg)*K*fs/fc;   %the number of samples for one diverse time period.
mm=[];
ff=[];
p_signal=0.5;
%snr=[-12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5];
%snr=[1 2 3 4 5 6 7 8];
snr=[-6 -5 -4 -3 -2 -1 0 1 2];
%Nblocks = [40 40 80 80 80 160 160 160 20 20 20 40 40 40 80 80 80 160 160 160 160 160 160];
Nblocks = [20 20 20 20 20 20 20 20 20 20 20 20 20 20 20];
[m n]=size(snr);
pe=[];
mppskfilter=[];

Nsymbols = time/20;   % the time period of baseband signal is 20ms
 b =[1  -1.414423081535762   1];     %fs/fc=8
 a =[1  -3.987007081485477   7.949450640987077  -9.385437352983406   7.009070118114206  -3.099682963193266   0.685515443313961];
 
 %% generate the MPPSK signals
   symbol_sample_stream_num=Nsymbols*n1;
        mppsk=zeros(Nsymbols,n1);
        s_tx =[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 randint(1,Nsymbols-85,[0 M - 1] )];  %[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 randint(1,Nsymbols-25,[0 M - 1])];
       i=1;
       snr(i)=0;   %-6    %6;
         fprintf('SNR = %2gdB\n',snr(i));
   %% calculate the noise and the SNR     
   %      snr_lin = 10^(snr(i)/10);
   % p_noise = p_signal/snr_lin;        
   % a_noise = sqrt(p_noise);
   % noise=a_noise*randn(1,symbol_sample_stream_num);
   % calculae the pulse noise and add it to the carrier
    SNRpbit=10.^(snr(i)/10);
    alpha=1.5;
 %   alpha=1.2;
%SNRpbit=10.^(dB/10);                % Eb/No conversion from dB to decimal
    R=1;% R=k/n1;                     % code rate 编码率
    No=1./SNRpbit;
    %由信噪比计算得到噪声产生的参数值
        S0=1/(sqrt(7.12*SNRpbit*R));        %3.56
        gamma=((1.78*S0)^alpha)/1.78;
        scale=gamma^(1/alpha);
       
        %调制和加噪 0->1,1->-1
      %  zb=-2*en_msg+1;   %编码这步怎么没有？encoder文件也打不开                      
      %  zb=(zb + stblrnd(alpha,0,scale,0,1,n));  %n个数编码之后为k个数   
 %   delta=sqrt(No/2)/sqrt(2);

    pulsenoise =0;   %stblrnd(alpha,0,scale,0,[1,symbol_sample_stream_num]);
   
    %% prepare the matrix and array for the signals
    t0 = (0:1:n1 - 1)/fs;
    f0 = A*sin(2*pi*fc*t0);    
    t1 = (0:1:n2 - 1)/fs;
    f1 = B*sin(2*pi*fc*t1 + thita);
    
    err = 0;
    s_rcv=[];
    s_tx_total=[];
    zf1 = zeros(max(length(b),length(a))-1,1);
    zf2 = zeros(max(length(b),length(a))-1,1);
%    zf_env1 = zeros(max(length(Num_low),1)-1,1);
%    zf_env2 = zeros(max(length(Num_low),1)-1,1);
%    zf_shape=zeros(max(length(b_shape),length(a_shape))-1,1);
         
          for num = 1:Nsymbols
            % the whole diverse do not need to seperate 0 and 1 2 3.... 
           % if( s_tx(num) == 0)       
 
           %     mppsk(num, :) = f0;                      
           % else                                        
                k = s_tx(num);                    
                head = f0(1:(k)*K*(fs/fc));     % head = f0(1:(k - 1)*K*(fs/fc));        
                middle=f1;
                tail = f0(1:(N - (k - rg+1)*K)*(fs/fc));    %tail = f0(1:(N - (k - rg)*K)*(fs/fc));   
                mppsk(num, :) = [head middle tail];     
            %end
          end
           mppsk = reshape(mppsk', 1, symbol_sample_stream_num);
     
         
        y_mppsk=mppsk+pulsenoise;       %+noise;    
        data = y_mppsk  ;  %reshape(y_mppsk', 1, symbol_sample_stream_num);
        
        
       

     
    %--- Initialization ---------------------------------------------------
    figure(100);
    
    
    clf(100);
    
    timeScale = 0 : 1/settings.samplingFreq : 40e-3;    
    
    %--- Time domain plot -------------------------------------------------
    subplot(2, 2, 1);
   % plot(1000 * timeScale(1:round(samplesPerCode/50)), ...
    %     data(1:round(samplesPerCode/50)));
      plot(1000* timeScale(1:round(samplesPerCode*40)), ...
         data(1:round(samplesPerCode*40)));
     
     
    axis tight;
    grid on;
    title ('Time domain plot');
    xlabel('Time (ms)'); ylabel('Amplitude');
    
    %--- Frequency domain plot --------------------------------------------
    subplot(2,2,2);
    pwelch(data-mean(data));   %, 16384, 1024, 2048, settings.samplingFreq/1e6)
    
    axis tight;
    grid on;
    title ('Frequency domain plot');
    xlabel('Frequency (MHz)'); ylabel('Magnitude');
    
    %--- Histogram --------------------------------------------------------
    subplot(2, 2, 3.5);
    hist(data, -128:128)
    
    dmax = max(abs(data)) + 1;
    axis tight;
    adata = axis;
    axis([-dmax dmax adata(3) adata(4)]);
    grid on;
    title ('Histogram'); 
    xlabel('Bin'); ylabel('Number in bin');

end % if (fid > 0)
