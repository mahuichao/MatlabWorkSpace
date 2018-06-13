function [FBFOutput,Output] = TF_GSC(x,fs,DesAng,NFFT,mu,forgottor)
% 函数说明：此版本为没有加入方向偏至F的版本
% 参数说明 mu,forgottor,threhod,fs,desAngle
% x:输入信号 M*L 默认是1000个数据进入一次
% win:窗函数
% overlap:步长
% NFFT:帧长
% Iden_NFFT:
% mu:梯度下降因子
% forgottor: 遗忘因子，Pest求取时使用
% fs:采样率
% DesAng:期望入射角  目前没有用到

%MArray = 0.0425*[1,0;1/2,sqrt(3/4);-1/2,sqrt(3)/2;-1,0;-1/2,-sqrt(3)/2;1/2,-sqrt(3)/2];
MArray=0.05*[-2.5;-1.5;-0.5;0.5;1.5;2.5];
% MArray=0.05*[-2.5,0;-1.5,0;-0.5,0;0.5,0;1.5,0;2.5,0];
%MArray=0.05*[0;1;2;3;4;5];
[MicNum,DataLength]=size(x);
window=hanning(NFFT);
overlap=NFFT/2;
N=fix((DataLength)/overlap)*overlap+overlap;
X=[x zeros(MicNum,N-DataLength)];  
nblocks=length(X)/overlap;
segs=13;
HThrehod=segs*NFFT; 
HSegments=fix(DataLength/HThrehod); 
c=340;
Output=zeros(N,1);
FBFOutput=zeros(N,1);
FBFOutbuff=zeros(NFFT/2,1);
H=zeros(MicNum,NFFT);
xbuf=zeros(MicNum,NFFT);
Xbuf=zeros(MicNum,NFFT);
lastOne=0;
W0=zeros(MicNum,NFFT);
YFBF=zeros(NFFT,1);
Y=zeros(NFFT,1);
HStd=zeros(MicNum,NFFT);
G=zeros(MicNum-1,NFFT);
outbuf=zeros(NFFT/2,1);
Pest=zeros(NFFT,1);
Freband=linspace(0,fs/2,NFFT/2+1);
timeDelay=zeros(MicNum,NFFT);
BM=zeros(MicNum,MicNum-1,NFFT);
%% steer vector
for freIdx=1:length(Freband)
    coefficient = 2 * pi * Freband(freIdx)/c;

     h=cosd(DesAng); % 模拟
%     h= [cosd(DesAng),sind(DesAng)].';
    tao = MArray*h;
    timeDelay(:,freIdx)=exp(j*coefficient*tao);
end


for m=1:MicNum
    tmpDelay=timeDelay(m,1:(NFFT/2+1)).';
    tmpConjD=timeDelay(m,2:NFFT/2).';
    timeDelay(m,:)=[tmpDelay;flipud(conj(tmpConjD))].';
end




%% 传统波束形成加权向量
% % w_con = zeros(MicNum,NFFT);
% % for FreIdx = 1 : length(Freband)
% %     coefficient = 2 * pi * Freband(FreIdx)/c;
% % %     h = [cosd(DesAng),sind(DesAng)].';
% %     h=cosd(DesAng); % 模拟
% %     tao = MArray*h;
% %     w_con(:,FreIdx) = exp(j*coefficient*tao)/MicNum;
% % end;
% % 
% % for m=1:MicNum
% %     tmpW_con=w_con(m,1:NFFT/2+1).';
% %     tmpW_conj=w_con(m,2:NFFT/2).';
% %     w_con(m,:)=[tmpW_con;flipud(conj(tmpW_conj))];
% % end
% 
% %% 自我测试
% TMH=Plot_H(X,NFFT,segs,window,timeDelay);
% n=1;
% subplot(331);
% 
% angle_01=atand(imag(TMH(1,(n-1)*NFFT+1:n*NFFT))./real(TMH(1,(n-1)*NFFT+1:n*NFFT)));
% % angle_01=atand(imag(TMH(1,1:end))./real(TMH(1,1:end)));
% plot(1:512,angle_01);
% % plot(1:length(TMH(1,:)),angle_01);
% title('Pic Mic1');
% subplot(332);
% angle_02=atand(imag(TMH(2,(n-1)*NFFT+1:n*NFFT))./real(TMH(2,(n-1)*NFFT+1:n*NFFT)));
% % angle_02=atand(imag(TMH(2,1:end))./real(TMH(2,1:end)));
% plot(1:512,angle_02);
% % plot(1:length(TMH(1,:)),angle_02);
% % polarplot(TMH(2 ,41),'*');
% title('Pic Mic2');
% subplot(333);
% angle_03=atand(imag(TMH(3,(n-1)*NFFT+1:n*NFFT))./real(TMH(3,(n-1)*NFFT+1:n*NFFT)));
% % angle_03=atand(imag(TMH(3,1:end))./real(TMH(3,1:end)));
% plot(1:512,angle_03);
% % plot(1:length(TMH(1,:)),angle_03);
% % polarplot(TMH(3 ,41),'*');
% title('Pic Mic3');
% subplot(334);
% angle_04=atand(imag(TMH(4,(n-1)*NFFT+1:n*NFFT))./real(TMH(4,(n-1)*NFFT+1:n*NFFT)));
% % angle_04=atand(imag(TMH(4,1:end))./real(TMH(4,1:end)));
% plot(1:512,angle_04);
% % plot(1:length(TMH(1,:)),angle_04);
% % polarplot(TMH(4 ,41),'*');
% title('Pic Mic4');
% subplot(335);
% angle_05=atand(imag(TMH(5,(n-1)*NFFT+1:n*NFFT))./real(TMH(5,(n-1)*NFFT+1:n*NFFT)));
% % angle_05=atand(imag(TMH(5,1:end))./real(TMH(5,1:end)));
% plot(1:512,angle_05);
% % plot(1:length(TMH(1,:)),angle_05);
% % polarplot(TMH(5 ,41),'*');
% title('Pic Mic5');
% subplot(336);
% angle_06=atand(imag(TMH(6,(n-1)*NFFT+1:n*NFFT))./real(TMH(6,(n-1)*NFFT+1:n*NFFT)));
% % angle_06=atand(imag(TMH(6,1:end))./real(TMH(6,1:end)));
% plot(1:512,angle_06);
% % plot(1:length(TMH(1,:)),angle_06);
% % polarplot(TMH(6 ,41),'*');
% title('Pic Mic6');


%% 算法开始   
for n=1:nblocks
  
    %% TF估计阶段
    % 查看目前在哪段H上，使用相应的H
    if(mod(n*overlap,HThrehod)~=0)
        HSIndex=fix(n*overlap/HThrehod)+1; 
    else
        HSIndex=fix(n*overlap/HThrehod);
    end
    
    if(lastOne~=HSIndex) % 并不是每次进来一帧 就需要运行一次
        if(HSIndex*HThrehod<=DataLength)  % 当超出索引的时候，我们就让H和W0等于上一次计算得到的值。
            % H=TF_Estimation(X(:,(HSIndex-1)*HThrehod+1:HSIndex*HThrehod),NFFT,13,window); 
            % aaaa=X(:,(HSIndex-1)*HThrehod+1:HSIndex*HThrehod);
            H=TF_Estimation_TD(X(:,(HSIndex-1)*HThrehod+1:HSIndex*HThrehod),NFFT,segs,window,timeDelay); 
        
%             for m=1:MicNum
%                 for freIdx=1:NFFT
%                     if H(m,freIdx)>1
%                         H(m,freIdx)=1;
%                     end
%                 end
%             end
%            FreBand=linspace(-fs/2,fs/2,NFFT);
%            plot(FreBand,H(2,:));
            
            
%             ifftH=ifft(H.');
%             ifftH=ifftH(1:181,:);
%             HStd=fft(ifftH,NFFT);
%             f=16000/2*linspace(-1,1,512);
%             plot(f,HStd(:,2));
            HStd=H.';
            lastOne=HSIndex;
           %% W0阶段
            for freIdx=1:NFFT
                W0(:,freIdx)=HStd(freIdx,:)./norm(HStd(freIdx,:)).^2    ; % MicNum*NFFT   .*timeDelay(:,freIdx).' 
            end
            
            
           %% Blocking Matrix
   
            for freIdx=1:NFFT
                for row=1:MicNum
                    for col=1:MicNum-1
                        if(row==1)
                            BM(row,col,freIdx)=-conj(HStd(freIdx,col+1));
                        elseif(row~=1 && row-1==col)
                            BM(row,col,freIdx)=1;
                        end
                    end
                end
            end
            
        end
       
    end
    
   %% GSC Beamforming阶段
    for m=1:MicNum
        xbuf(m,NFFT/2+1:end)=X(m,(n-1)*overlap+1:n*overlap);
        Xbuf(m,:)=fft(xbuf(m,:)'.*window) ; % .*rectwin(NFFT)
        Xbuf(m,:)=Xbuf(m,:).*timeDelay(m,:) ;
    end

    %% YFBF  求解
    for freIdx=1:NFFT
        YFBF(freIdx,1)=W0(:,freIdx)'*Xbuf(:,freIdx); % NFFT*MicNum  *  MicNum*1=NFFT*1
    end
    

    
    
   % YFBF(:,1)=Y_con(:,1);
   % YFBF(:,1)=Xbuf(2,:).';
    %===============测试开始===================
%       YBF_time=real(ifft(YFBF));
%       FBFOutput((n-1)*NFFT/2+1:n*NFFT/2)=FBFOutbuff+YBF_time(1:NFFT/2); %
%       FBFOutbuff=YBF_time(NFFT/2+1:end);
%       xbuf(:,1:NFFT/2)=xbuf(:,NFFT/2+1:end);
    %===============测试结束===================
%      continue;
     

    
    U=zeros(MicNum-1,NFFT);
    for freIdx=1:NFFT
        U(:,freIdx)=BM(:,:,freIdx)'*Xbuf(:,freIdx); % (M-1) *M  *M*1=(M-1)*1 
    end
    %% NC
    if n==1
        Y=YFBF; % G的初始值设为0带入论文 （21）和（22） 之间的公式
    else
        for freIdx=1:NFFT
            Y(freIdx,1)=YFBF(freIdx,1)-G(:,freIdx)'*U(:,freIdx); % 1*1 - 1*(M-1)*(M-1)*1=1
        end
    end
    %% LMS
 
    for freIdx=1:NFFT
        Pest(freIdx,1)=forgottor*Pest(freIdx,1)+(1-forgottor)*norm(Xbuf(2:end,freIdx)).^2;
        G(:,freIdx)=G(:,freIdx)+mu*(U(:,freIdx)*conj(Y(freIdx,1)))/Pest(freIdx,1); % (M-1)*1  G在不断的增大
    end

%     for m=1:MicNum-1
%         for freIdx=1:NFFT
%             if G(m,freIdx)>3e-6;
%                G(m,freIdx)=3e-6;
%             end
%         end
%     end
    
  %  Freband=linspace(-fs/2,fs/2,NFFT);
  % plot(Freband,G(2,:));
   G_time=ifft(G.')';
   G_time=G_time(:,1:251); % 对应论文中的251个系数的FIR滤波器(不过效果并不是很大)
   G=fft(G_time.',NFFT).';     
    
    %% 通用计算部分
    
    Y_time=real(ifft(Y));
    Output((n-1)*NFFT/2+1:n*NFFT/2)=outbuf+Y_time(1:NFFT/2); %
    outbuf=Y_time(NFFT/2+1:end);
    xbuf(:,1:NFFT/2)=xbuf(:,NFFT/2+1:end);
end
%Output=Output.*4; % 由于加窗带来的幅度影响，所以我们*2来消除hanning窗的0.5带来的影响(幅度相等原则)
%name=['origin_last_02_' num2str(mu) '_' num2str(forgottor) '.wav']
%audiowrite(name,Output(:,1),fs);
end