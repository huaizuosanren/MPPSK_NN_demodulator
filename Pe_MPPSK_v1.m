clc;
clear all;
close all;

M=128;
K = 2;
N =260;
A=1;
B=1;
rg=0;
thita=pi;


fc=981e3;     
fs=fc*20; 
n1=N*fs/fc;
n2=(1-rg)*K*fs/fc;
mm=[];
ff=[];
p_signal=0.5;
%snr=[-12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5];
%snr=[1 2 3 4 5 6 7 8];
snr=[ -2 -1 0 1 2]
Nblocks=[20 20 20 20 20 20 20 20 20 20 20 20 20 20 20]
[m n]=size(snr);
pe=[];
mppskfilter=[];

Nsymbols =5000;

if fs/fc==10
     b=[1,-1.618092409933249,0.99990000250000044];
     a=[1,-4.5620074920961651,9.5862839416819483,-11.566980661101638,8.4523528839743243,-3.5467147693005732,0.6855154433139603];
     Hd=lowfilter_eb_10;
     b_env=conv( conv( Hd.sosMatrix(1,1:3),Hd.sosMatrix(2,1:3) ),Hd.sosMatrix(3,1:3) )*Hd.ScaleValues(1)*Hd.ScaleValues(2)*Hd.ScaleValues(3)*Hd.ScaleValues(4);
     a_env=conv( conv( Hd.sosMatrix(1,4:6),Hd.sosMatrix(2,4:6) ) ,Hd.sosMatrix(3,4:6));
 elseif fs/fc==8
     b =[1  -1.414423081535762   1];
     a =[1  -3.987007081485477   7.949450640987077  -9.385437352983406   7.009070118114206  -3.099682963193266   0.685515443313961];
     Hd=lowfilter_eb_8;
     b_env=conv( conv( Hd.sosMatrix(1,1:3),Hd.sosMatrix(2,1:3) ),Hd.sosMatrix(3,1:3) )*Hd.ScaleValues(1)*Hd.ScaleValues(2)*Hd.ScaleValues(3)*Hd.ScaleValues(4);
     a_env=conv( conv( Hd.sosMatrix(1,4:6),Hd.sosMatrix(2,4:6) ) ,Hd.sosMatrix(3,4:6));
 elseif fs/fc==6
     b =[1  -1.000272874264982   1];
     a =[1  -2.813475240872029   5.282770438671800  -5.776182128097824   4.640779780443355  -2.171483139266619   0.678084436470381];
     Hd=lowfilter_eb_6;
     b_env=conv( conv( Hd.sosMatrix(1,1:3),Hd.sosMatrix(2,1:3) ),Hd.sosMatrix(3,1:3) )*Hd.ScaleValues(1)*Hd.ScaleValues(2)*Hd.ScaleValues(3)*Hd.ScaleValues(4);
     a_env=conv( conv( Hd.sosMatrix(1,4:6),Hd.sosMatrix(2,4:6) ) ,Hd.sosMatrix(3,4:6));
 elseif fs/fc==4
     b=[1,-0.00039265358727080812,1];
     a=[1,0.001712686483139045 ,2.6563009792332659 ,0.0030264651110280486 ,2.3418192253785417,0.0013315213654196743 ,0.68551544331395908];
     Hd=lowfilter_eb_4; 
     b_env=conv(Hd.sosMatrix(1,1:3),Hd.sosMatrix(2,1:3))*Hd.ScaleValues(1)*Hd.ScaleValues(2)*Hd.ScaleValues(3);
     a_env=conv(Hd.sosMatrix(1,4:6),Hd.sosMatrix(2,4:6));
elseif fs/fc==100
     b = [1,-1.9960599659357392, 1];
     a = [1,-3.8715335529480597 ,  5.6264550392653057 ,-3.6377858431552683,  0.8828932613760001];
     b_env=[0.00000000000000019979744992511669,0.0000000000000017981770493260502 ,0.0000000000000071927081973042008 ,0.000000000000016782985793709801  ,0.0000000000000251744786905647 ,0.0000000000000251744786905647 ,0.000000000000016782985793709801 ,0.0000000000000071927081973042008 ,0.0000000000000017981770493260502,0.00000000000000019979744992511669];
     a_env=[1,  -8.7901082724397721,  34.342845035539966,  -78.275384915157446, 114.69870007631428,-112.05472119025141,  72.986069384142766, -30.562822163384268,   7.4660809754742745,  -0.81065893023827074];
     load mppsk_in100;
     load s_tx_f100;
     num_round=21;
     fprintf('经过成形滤波  采样率：100倍\n \n');
elseif fs/fc==60
     b = [1,-1.9890562794070863 , 1];
     a = [1,-3.7785625564254337 , 5.3689524801538155 ,-3.3996860888926554 , 0.80951407289999999];
     b_env=[0.000000000000018736618234530585,0.00000000000016862956411077527,0.00000000000067451825644310117,0.0000000000015738759317005695,0.0000000000023608138975508541,0.0000000000023608138975508537,0.0000000000015738759317005691,0.00000000000067451825644310117,0.00000000000016862956411077527,0.000000000000018736618234530585];
     a_env=[1,-8.6497187354147815,33.258873379430888,-74.613054426804382,107.62679640162835,-103.51838880608199,66.390362131290445,-27.377117501346461,6.5866846140872859,-0.70443705677975388];
elseif fs/fc==20
    b=[1,-1.9021496572560159 ,1];
    a=[1,-3.6512241163814698 , 5.1727286626648894 ,-3.3577825365961242, 0.84572301542400019];
%     b_env=[0.00000000030617715438660676  0.0000000027555943894794609  0.000000011022377557917844  0.000000025718880968474968  0.000000038578321452712451  0.000000038578321452712451  0.000000025718880968474968  0.000000011022377557917842  0.0000000027555943894794609  0.00000000030617715438660676];
%     a_env=[1  -7.9314941219361126   28.016432168224561  -57.840833562173017   76.909993309948646  -68.299810046586828   40.505579522030544  -15.468452202861751    3.4513742035511221   -0.34278911343445534];
%     b_shape=[0.132499936880248  -0.000000000000888  -0.264999873759386  -0.000000000002665   0.132499936881136]*1e-3;
%     a_shape=[1.000000000000000  -3.773188120335216   5.526651408008984  -3.711763047740910   0.967708444925353];
%     SOS=[1,0,-1,1, -1.8836900259519243, 0.98754643299250611;
%         1,0,-1,1,-1.8973481693610572 , 0.98832263073963533;
%         1,0,-1,1,-1.8793863998698668 , 0.97601099599076246;
%         ];
    load SOS;
    load num_low;
    [b_shape,a_shape]=sos2tf(SOS,G);
    fprintf('未经滤波  采样率：20倍\n \n');
    

elseif fs/fc==30
    b=[1,-1.9563199437524734 ,1];
    a=[1,-3.7943813789556295, 5.4787341498319906,-3.5652918088198042, 0.88289326137599999];
    b_env=[0.00000000052344979143157986, 0.0000000041875983314526389, 0.000000014656594160084236 ,0.000000029313188320168475 ,0.00000003664148540021059 ,0.000000029313188320168475 ,0.000000014656594160084236 ,0.0000000041875983314526389 ,0.00000000052344979143157986];
    a_env=[1, -7.2586867437728975, 23.08286154866521,-42.000652224794585, 47.824926893700308,-34.895140409225043 , 15.931981890687625 , -4.1613458133859744 ,  0.47605499212850316];
    load mppsk_in30;
    load s_tx_f30;
    fprintf('未经滤波  采样率：30倍\n \n');
    num_round=61;
    
elseif fs/fc==40
    b=[1,-1.9753890794728897 ,1];
    a=[1,-3.7128942681928949 , 5.2060421878608025 ,-3.2660400172273425 , 0.77378060390399828];
    b_env=[0.0000000000006773692970349954,0.0000000000060963236733149589,0.000000000024385294693259836 ,0.000000000056899020950939623 ,0.000000000085348531426409441 ,0.000000000085348531426409441 ,0.000000000056899020950939629 ,0.000000000024385294693259839 ,0.0000000000060963236733149589,0.0000000000006773692970349954];
    a_env=[1, -8.4732195544588382, 31.923740256506051,-70.192581084425484, 99.260075827315717,-93.616949641020938, 58.888064699180219,-23.822935343308693,  5.6241571866416793, -0.59035234608287479];
    load mppsk_in40;
    load s_tx_f40;
    fprintf('未经滤波  采样率：40倍\n \n');
    num_round=41;

elseif fs/fc==50
    b=[1,-1.9842403563983937,1];
    a=[1,-3.7791412146145271, 5.3798937602070085,-3.4187547501064568, 0.81836991104400025];
    b_env=[0.000000000000094215049270986409,0.00000000000084793544343887782,0.0000000000033917417737555113,0.0000000000079140641387628583,0.000000000011871096208144287,0.000000000011871096208144287,0.0000000000079140641387628583,0.0000000000033917417737555105,0.00000000000084793544343887782,0.000000000000094215049270986409];
    a_env=[1, -8.5792806236814734, 32.722357591877426,-72.824657074551297,104.21931137343614,-99.45948947962188, 63.295189286596404,-25.901540399133847,  6.1845982874205747, -0.65648896229381382];
    load mppsk_in50;
    load s_tx_f50;
    fprintf('未经滤波  采样率：50倍\n \n');
    num_round=41;
    
end

for(i=1:n)
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
        mppsk=zeros(Nsymbols,n1);    %n1 is the sample numbers for every baseband signals
        s_tx = randint(1,Nsymbols,[0 M - 1]);
%         if j==1
%             for jjj=1000;Nsymbols
%                 s_tx[jjj]=mod(jjj,M);
%             end
%         end
        %得到调制信号并加噪
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
        
        [mppsk,zf_shape]=filter(b_shape,a_shape,mppsk);
%         
        noise=a_noise*randn(1,symbol_sample_stream_num);
        y_mppsk=mppsk+noise;
        
        timeScale = 0 : 1/fs : 0.05; 
        
    %    figure (100)
    %    subplot(2,1,1);
    %    plot(timeScale,mppsk(1:length(timeScale)));
    %    subplot (2,1,2);
    %    plot(timeScale,y_mppsk(1:length(timeScale)));
    %    title('sending signals');
        
        %冲击滤波
        [mppskfilter,zf1]=filter(b,a,mppsk,zf1);
        [y_mppskfilter,zf2]=filter(b,a,y_mppsk,zf2);
   %     figure (101)
   %     subplot(2,1,1);
   %     plot(timeScale,mppskfilter(1:length(timeScale)));
   %     subplot (2,1,2);
   %     plot(timeScale,y_mppskfilter(1:length(timeScale)));
   %     title('pulsefiltered signals');
        %******************取绝对值******************************
        mppskfilter=abs(mppskfilter);
        y_mppskfilter=abs(y_mppskfilter);
        mppsk_length=length(mppskfilter);
    %    figure (102)
    %    subplot(2,1,1);
    %    plot(timeScale,mppskfilter(1:length(timeScale)));
    %    subplot (2,1,2);
    %    plot(timeScale,y_mppskfilter(1:length(timeScale)));
    %    title('abs signals');
% %         %*****************取包络***********************************
        [mppsk_baoluo,zf_env1]=filter(Num_low,1,mppskfilter,zf_env1);
        [y_mppsk_baoluo,zf_env2]=filter(Num_low,1,y_mppskfilter,zf_env2);
    %    figure (103)
    %    subplot(2,1,1);
    %    plot(timeScale,mppsk_baoluo(1:length(timeScale)));
    %    subplot (2,1,2);
    %    plot(timeScale,y_mppsk_baoluo(1:length(timeScale)));
    %    title('baoluo signals');

%         ******得到发送序列中各个值的个数，以便于确定阈值*********


         if j==1
            ind0=find(s_tx(1001:Nsymbols)==0);                                      
            Nsymbols0=length(ind0);
            mppsk_baoluo0=zeros(n1,Nsymbols0);          
            y_mppsk_baoluo0=zeros(n1,Nsymbols0); 
            for l=1:Nsymbols0
                mppsk_baoluo0(:,l)=mppsk_baoluo((1001+ind0(l)-2)*n1+1:(1001+ind0(l)-1)*n1)';
                y_mppsk_baoluo0(:,l)=y_mppsk_baoluo((1001+ind0(l)-2)*n1+1:(1001+ind0(l)-1)*n1)';
            end
            for l=1:M-1
                ind1=find(s_tx(1001:Nsymbols)==l);
                Nsymbols1=length(ind1);
                mppsk_baoluo1=zeros(n1,Nsymbols1);
                for il=1:Nsymbols1
                    mppsk_baoluo1(:,il)=mppsk_baoluo((1001+ind1(il)-2)*n1+1:(1001+ind1(il)-1)*n1)';
                    y_mppsk_baoluo1(:,il)=y_mppsk_baoluo((1001+ind1(il)-2)*n1+1:(1001+ind1(il)-1)*n1)';
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
          %  Vd_jf_noise
         else 
              s_rcv=zeros(1,Nsymbols);
              for t=1:Nsymbols
                  s_rcv_al=0;
                  s_rcv_tmp=1;           
                  symbol_tmp=y_mppsk_baoluo((t-1)*n1+1:t*n1);            
                  for l=1:M-1
                      y_mppsk_baoluo_jf1=sum(symbol_tmp(det_k(l)-3:det_k(l)+3));
                      if s_rcv_al < y_mppsk_baoluo_jf1
                          s_rcv_al = y_mppsk_baoluo_jf1;
                          s_rcv_tmp = l;
                      end 
                  end
                  if (s_rcv_al > Vd_jf_noise(s_rcv_tmp))
                      s_rcv(t) = s_rcv_tmp;
                  else
                      s_rcv(t)=0;
                  end
              end
             
              err=err+sum(s_tx~=s_rcv);
              
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
    pe(i)=err/(Nblocks(i)*Nsymbols)
    num0
    num1
    num2
    num3
end

semilogy(snr,pe);
grid on;
str= sprintf('pe%d',fs/fc);
save(str,'pe');
str=sprintf('snr%d',fs/fc);
save(str,'snr');