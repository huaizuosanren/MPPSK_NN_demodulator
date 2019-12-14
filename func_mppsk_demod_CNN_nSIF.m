function [snr,pe]= func_mppsk_demod_CNN_nSIF(M,N,K,code)


%M= 8  %2  %4  %2  %16   %64    %128;     %128进制
%K = 57 %48   % 20  %2  %99  %198  %24    %6     %3; %2;    % every diverse is 2T
%N = 399 %400  %180  %10  %400;   %260;   % every base band signal is 260T
M
N
K
code
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
%snr=[-5 -10 -8 -6 -4 -2 0];
%snr=[-2 0 2 4 6 8];
snr=[-6 -4 -2 0 2 4 6];
%snr=[-6 -5 -4 -3 -2 -1 0 1 2];
%snr=[2 4 6 8 10 12 14];
Nblocks = [40 40 80 80 80 160 160 160 160 160 160 160 160];
%$Nblocks = [20 20 20 20 20 20 20 20 20 20 20 20 20 20 20];
[m n]=size(snr);
pe=[];
mppskfilter=[];

Nsymbols = 5000;  %20000;  %5000;

 % b=[1,-1.9021496572560159 ,1];     % fs/fc=20
 %   a=[1,-3.6512241163814698 , 5.1727286626648894 ,-3.3577825365961242, 0.84572301542400019];
 
   b =[1  -1.414423081535762   1];     %fs/fc=8
     a =[1  -3.987007081485477   7.949450640987077  -9.385437352983406   7.009070118114206  -3.099682963193266   0.685515443313961];
    % Hd=lowfilter_eb_8;
   %  b_env=conv( conv( Hd.sosMatrix(1,1:3),Hd.sosMatrix(2,1:3) ),Hd.sosMatrix(3,1:3) )*Hd.ScaleValues(1)*Hd.ScaleValues(2)*Hd.ScaleValues(3)*Hd.ScaleValues(4);
   %  a_env=conv( conv( Hd.sosMatrix(1,4:6),Hd.sosMatrix(2,4:6) ) ,Hd.sosMatrix(3,4:6));
   
    Hd=tianxian; % the tianxian filter
 
 
      load SOS;
    load num_low;
    [b_shape,a_shape]=sos2tf(SOS,G);
    fprintf('未经滤波  采样率：8倍\n \n');
    
    
%% initialize the variables---------------------------------------------------------------------
        symbol_sample_stream_num=Nsymbols*n1;
        mppsk=zeros(Nsymbols,n1);
        s_tx =[];   %randint(1,Nsymbols,[0 M - 1]);
        
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
   %**********calculate the amplifier of the gaosi noise*****
   % snr_lin = 10^(snr(i)/10);
   % p_noise = p_signal/snr_lin;        
   % a_noise = sqrt(p_noise);
    
     %**************calculate the puse noise************* 
    SNRpbit=10.^(snr(i)/10);
    No=0.5./SNRpbit;
    delta=sqrt(No)/sqrt(2);
    alpha=1.5;
    pulsenoise = stblrnd(alpha,0,delta,0,[1,symbol_sample_stream_num]);
    
    t0 = (0:1:n1 - 1)/fs;
    f0 = A*sin(2*pi*fc*t0);    
    t1 = (0:1:n2 - 1)/fs;
    f1 = B*sin(2*pi*fc*t1 + thita);
    
    err = 0;
    s_rcv=[];
    s_tx_total=[];
    zf1 = zeros(max(length(b),length(a))-1,1);
    zf2 = zeros(max(length(b),length(a))-1,1);
    zf_env1 = zeros(max(length(Num_low),1)-1,1);
    zf_env2 = zeros(max(length(Num_low),1)-1,1);
    zf_shape=zeros(max(length(b_shape),length(a_shape))-1,1);
    
    
    
     for j = 1:1+Nblocks(i)
        if(mod(j-1,20)==19)
            fprintf('需循环次数:%d 完成：%f \n',Nblocks(i),(j-1)/Nblocks(i));
        end
        
        symbol_sample_stream_num=Nsymbols*n1;
        mppsk=zeros(Nsymbols,n1);
        s_tx = randi([0 M - 1],1,Nsymbols);   %randint(1,Nsymbols,[0 M - 1]);
    
 %% 得到调制信号并加噪
        for num = 1:Nsymbols
           %  the whole diverse do not need to seperate 0 and 1 2 3.... 
            if( s_tx(num) == 0)       
 
                mppsk(num, :) = f0;                      
            else                                        
                k = s_tx(num);                    
                head = f0(1:(k)*K*(fs/fc));     % head = f0(1:(k - 1)*K*(fs/fc));        
                middle=f1;
                tail = f0(1:(N - (k - rg+1)*K)*(fs/fc));    %tail = f0(1:(N - (k - rg)*K)*(fs/fc));   
                mppsk(num, :) = [head middle tail];     
            end
        end
        mppsk = reshape(mppsk', 1, symbol_sample_stream_num);
        
        %*******************Pass Through The tianxian***************************
          mppsk=filter(Hd,mppsk);
        
   %     [mppsk,zf_shape]=filter(b_shape,a_shape,mppsk);
%         
     %******************calcaulate the gaosi
     %noise***************************
     %   noise=a_noise*randn(1,symbol_sample_stream_num);
     
     %   y_mppsk=mppsk;   %+pulsenoise;    %+noise;    
         y_mppsk=mppsk+pulsenoise; 
     
     %----------限幅滤波--------------
     for mm=1:length (y_mppsk)
           if y_mppsk(mm)>=1.5
               y_mppsk(mm)=1;
           elseif y_mppsk(mm)<=-1.5
                 y_mppsk(mm)=-1;
           end
     end
     
  
 %        [mppskfilter,zf1]=filter(b,a,mppsk,zf1);
 %       [y_mppskfilter,zf2]=filter(b,a,y_mppsk,zf2);
        y_mppsk_baoluo=y_mppsk;     %=y_mppskfilter;
        
          %******************取绝对值******************************
      %  mppskfilter=abs(mppskfilter);
      %  y_mppskfilter=abs(y_mppskfilter);
      %  mppsk_length=length(mppskfilter);

% %         %*****************取包络***********************************
      %  [mppsk_baoluo,zf_env1]=filter(Num_low,1,mppskfilter,zf_env1);
      %  [y_mppsk_baoluo,zf_env2]=filter(Num_low,1,y_mppskfilter,zf_env2);

 %    end
%% ***********learn the signals' input and output under this SNR**************     
 % train_t=1:n1;
     if j==1   % && i==1
      net=patternnet(code);    %(10);
      net.trainParam.epochs=200;  %100;  % iterations
      net.trainParam.lr=0.001;  % learning ratio
      net.trainParam.goal=0.00004;  % mu biao wu cha
      net.trainParam.max_fail=8;
    %  net.trainParam.showWindow=false;
      y_mppsk_baoluo_tmp=reshape(y_mppsk_baoluo,n1,Nsymbols);
    %  y_mppsk_baoluo_tmp=(y_mppsk_baoluo_tmp)';
     % for t=1:Nsymbols
     %   y_mppsk_baoluo_tmp(1,:)=(y_mppsk_baoluo((t-1)*n1+1:t*n1))';
     %   y_mppsk_baoluo_tmp(2,:)=train_t';
     %   s_tx_tmp=zeros(Nsymbols,n1);
     % for t=1:Nsymbols
       % for s_tx_i=1:n1
     %     s_tx_tmp(t,:)=s_tx(t);
       % end
     % end
     
     s_tx_tmp=ind2vec(s_tx+1);
     s_tx_tmp=full(s_tx_tmp);
        net=train(net,y_mppsk_baoluo_tmp,s_tx_tmp);
     % end
 %     view net
      else
      y_mppsk_baoluo_tmp=reshape(y_mppsk_baoluo,n1,Nsymbols);
 %    y_mppsk_baoluo_tmp=(y_mppsk_baoluo_tmp)';
     % s_rcvtmp=zeros(Nsymbols,n1);
        %for t=1:Nsymbols
         % y_mppsk_baoluo_tmp(1,:)=(y_mppsk_baoluo((t-1)*n1+1:t*n1))';
         % y_mppsk_baoluo_tmp(2,:)=train_t';
          s_rcvtmp=net(y_mppsk_baoluo_tmp);  % s_rcvtmp is Nsymbols*n1
         % s_rcvtmp=(s_rcvtmp)';  % s_rcvtmp is n1*Nsymbols
         % rcv=(sum(s_rcvtmp)/3200);  %rcv is 1*Nsymbols
         rcv=vec2ind(s_rcvtmp);
          s_rcv= rcv-1;    %vec2ind(s_rcvtmp);
        %end
         err=err+sum(s_tx~=s_rcv);
                           for ii=1:Nsymbols
                  if s_tx(ii)~=s_rcv(ii)
                      if s_tx(ii)==1;
                          num1=num1+1;
%                           err_one(num1)=s_rcv(ii);
                      elseif s_tx(ii)==0
                          num0=num0+1;
                      elseif s_tx(ii)==2
                          num2=num2+1;
%                           err_two(num2)=s_rcv(ii);
                      elseif s_tx(ii)==3
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
         
   pe(i)=err/(Nblocks(i)*Nsymbols)
    num0
    num1
    num2
    num3
    clear y_mppsk_baoluo y_mppsk_baoluo_tmp y_mppsk;
   end
   
  %figure 1;
  %subplot(2,1,1);
  %plot(s_tx);
  %title('send out');
  %subplot(2,1,2);
  %plot(s_rcv);
  %title('receive')
  

  