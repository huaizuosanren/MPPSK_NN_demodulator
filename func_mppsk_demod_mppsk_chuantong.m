function [snr,pepulse]= func_mppsk_demod_mppsk_chuantong (snr,M,K,N,Vd,alpha)

%M= 8   %64  %2 %2   %3 %2  %64    %128;     %128jizhi
%K =57 % 48  %6  %8  %8  %133  %198  %6     %3; %2;    % every diverse is 2T
%N = 399 %400  %20   %400;   %260;   % every base band signal is 260T
A=1;     % amplifier is 1
B=1;
rg=0;   
thita=pi;  % the diverse angle is pi
M
K
N
Vd
alpha
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
[m n]=size(snr);
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
        err=zeros(1,n);
        
   for(i=1:n)
       % i=1;
     %   snr(i)=6;   %-6    %6;
         fprintf('SNR = %2gdB\n',snr(i));
    num0=0;
    num1=0;
    num2=0;
    num3=0;
    num4=0;
    err_three=zeros(1,500);
              err_one=zeros(1,300);
              err_two=zeros(1,500);
              err_four=zeros(1,500);
    snr_lin = 10^(snr(i)/10);
    p_noise = p_signal/snr_lin;        
    a_noise = sqrt(p_noise);
    
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
   %**************calculate the puse noise******************** 
    SNRpbit=10.^(snr(i)/10);
    No=0.5./SNRpbit;
    delta=sqrt(No)/sqrt(2);
%    alpha=1.5;
    
    
     for j = 1:1+Nblocks(i)
        if(mod(j-1,20)==19)
            fprintf('time to loop:%d complished£º%f \n',Nblocks(i),(j-1)/Nblocks(i));
        end
        
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
        mppsk=filter(Hd,mppsk);
        
%********caculate the pulse noise*************
    pulsenoise = stblrnd(alpha,0,delta,Vd,[1,symbol_sample_stream_num]);
%         
       % noise=a_noise*randn(1,symbol_sample_stream_num);
        y_mppsk=mppsk+pulsenoise;    
    %    mppsk=y_mppsk;
        
  %****************xianfulvbo****************************
       for mm=1:length (y_mppsk)
           if y_mppsk(mm)>=1.5+Vd
               y_mppsk(mm)=1+Vd;
           elseif y_mppsk(mm)<=-1.5+Vd
                 y_mppsk(mm)=-1+Vd;
           end
       end
    
  
         [mppskfilter,zf1]=filter(b,a,mppsk,zf1);
        [y_mppskfilter,zf2]=filter(b,a,y_mppsk,zf2);
        
          %******************get the abs******************************
        mppskfilter=abs(mppskfilter);
        y_mppskfilter=abs(y_mppskfilter);
        mppsk_length=length(mppskfilter);

% %         %*****************get the baoluo***********************************
        [mppsk_baoluo,zf_env1]=filter(Num_low,1,mppskfilter,zf_env1);
        [y_mppsk_baoluo,zf_env2]=filter(Num_low,1,y_mppskfilter,zf_env2);
%         ******get the number of every signals of send series, to define the yuzhi *********


         if j==1   % && i==1
            ind0=find(s_tx(((Nsymbols/5)+1):Nsymbols)==0);                     %find(s_tx(1001:Nsymbols)==0);                                      
            Nsymbols0=length(ind0);
            mppsk_baoluo0=zeros(n1,Nsymbols0);          
            y_mppsk_baoluo0=zeros(n1,Nsymbols0); 
            for l=1:Nsymbols0
                mppsk_baoluo0(:,l)= mppsk_baoluo((((Nsymbols/5)+1)+ind0(l)-2)*n1+1:(((Nsymbols/5)+1)+ind0(l)-1)*n1)';                    %mppsk_baoluo((1001+ind0(l)-2)*n1+1:(1001+ind0(l)-1)*n1)';
                y_mppsk_baoluo0(:,l)=y_mppsk_baoluo((((Nsymbols/5)+1)+ind0(l)-2)*n1+1:(((Nsymbols/5)+1)+ind0(l)-1)*n1)';                                             %y_mppsk_baoluo((1001+ind0(l)-2)*n1+1:(1001+ind0(l)-1)*n1)';
            end
            for l=1:M-1
                ind1=find(s_tx(((Nsymbols/5)+1):Nsymbols)==l);   %find(s_tx(1001:Nsymbols)==l);
                Nsymbols1=length(ind1);
                mppsk_baoluo1=zeros(n1,Nsymbols1);
                for il=1:Nsymbols1
                    mppsk_baoluo1(:,il)=mppsk_baoluo((((Nsymbols/5)+1)+ind1(il)-2)*n1+1:(((Nsymbols/5)+1)+ind1(il)-1)*n1)';        %mppsk_baoluo((1001+ind1(il)-2)*n1+1:(1001+ind1(il)-1)*n1)';
                    y_mppsk_baoluo1(:,il)=y_mppsk_baoluo((((Nsymbols/5)+1)+ind1(il)-2)*n1+1:(((Nsymbols/5)+1)+ind1(il)-1)*n1)';       %y_mppsk_baoluo((1001+ind1(il)-2)*n1+1:(1001+ind1(il)-1)*n1)';
                end
                det_k(l)=find(mppsk_baoluo1(:,1)==max(mppsk_baoluo1(:,1)));
                y_mppsk_baoluo1=reshape(y_mppsk_baoluo1,n1,[]);
     
                y_mppsk_baoluo_jf1=sum(y_mppsk_baoluo1(det_k(l)-3:det_k(l)+3,1:Nsymbols1));
                y_mppsk_baoluo_jf01=sum(y_mppsk_baoluo0(det_k(l)-3:det_k(l)+3,1:Nsymbols0));  
                signal_ave=sum(y_mppsk_baoluo_jf1(:))/Nsymbols1;
                noise_ave=sum(y_mppsk_baoluo_jf01(:))/Nsymbols0;
                Vd_jf_noise(l)=(signal_ave+noise_ave)/2*0.8;
                clear mppsk_baoluo1 y_mppsk_baoluo1 y_mppsk_baoluo_jf1 y_mppsk_baoluojf01
            end
           % det_k
         else 
              s_rcv=zeros(1,Nsymbols);
              for t=1:Nsymbols
                  s_rcv_al=0;     %0;
                  s_rcv_tmp=0;    %1;           
                  symbol_tmp=y_mppsk_baoluo((t-1)*n1+1:t*n1);            
                  for l=1:M-1
                      y_mppsk_baoluo_jf1=sum(symbol_tmp(det_k(l)-3:det_k(l)+3));
                      if s_rcv_al < y_mppsk_baoluo_jf1
                          s_rcv_al = y_mppsk_baoluo_jf1;
                          s_rcv_tmp = l;
                       %   s_rcv(t) = s_rcv_tmp;
                      end 
                  end
                  if (s_rcv_al > Vd_jf_noise(s_rcv_tmp))
                      s_rcv(t) = s_rcv_tmp;
                  else
                      s_rcv(t)=0;
                  end
              end
             
              err(i)=err(i)+sum(s_tx~=s_rcv);
              
              for ii=1:Nsymbols
                  if s_tx(ii)~=s_rcv(ii)
                      if s_tx(ii)==1;
                          num1=num1+1;
%                           err_one(num1)=s_rcv(ii);
                      elseif s_tx(ii)==0
                          num0=num0+1;
                      elseif s_tx(ii)==5
                          num2=num2+1;
%                           err_two(num2)=s_rcv(ii);
                      elseif s_tx(ii)==10
                          num3=num3+1;
%                           err_three(num3)=s_rcv(ii);
                      elseif s_tx(ii)==100
                          num4=num4+1;
%                           err_four(num4)=s_rcv(ii);
                      end
                  end
              end
         end
     end
         
   pepulse(i)=err(i)/(Nblocks(i)*Nsymbols)
   err
    num0
    num1
    num2
    num3
   end
 %% *****************plot the fiure*************************** 
  %--- Initialization ---------------------------------------------------
    figure(100);
    
    
    clf(100);
    
    timeScale = 0 : 1/fs : 0.01;    
    
    %--- Time domain plot -------------------------------------------------
    subplot(2, 2, 1);
   % plot(1000 * timeScale(1:round(samplesPerCode/50)), ...
    %     data(1:round(samplesPerCode/50)));
      plot(1000* timeScale, ...
         mppsk(1:length(timeScale)));
     
     
    axis tight;
    grid on;
    title ('Time domain plot');
    xlabel('Time (ms)'); ylabel('Amplitude');
    
    %--- Frequency domain plot --------------------------------------------
    subplot(2,2,2);
    pwelch(mppsk-mean(mppsk));   %, 16384, 1024, 2048, settings.samplingFreq/1e6)
    
    axis tight;
    grid on;
    title ('Frequency domain plot');
    xlabel('Frequency (MHz)'); ylabel('Magnitude');
    
    %--- Histogram --------------------------------------------------------
    subplot(2, 2, 3.5);
    hist(mppsk, -128:128)
    
    dmax = max(abs(mppsk)) + 1;
    axis tight;
    adata = axis;
    axis([-dmax dmax adata(3) adata(4)]);
    grid on;
    title ('Histogram'); 
    xlabel('Bin'); ylabel('Number in bin');
   

  %figure 1;
  %subplot(2,1,1);
  %plot(s_tx);
  %title('send out');
  %subplot(2,1,2);
  %plot(s_rcv);
  %title('receive')
 