close all; clc; clear all;

%% 首先我们模拟一下麦克风阵列，不如就设置为6个麦克风，间距为5cm，ULA也就是说成线性排列（不会的去看我csdn之前写的一片博客，有讲解模拟ULA的）
%% 6个麦克风，5cm。
microphone = ...
    phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 20e3]);

Nele = 6;
ula = phased.ULA(Nele,0.05,'Element',microphone);
c = 340;                 

%假设我们认为除了纯净语音外，还有人的笑声
ang_cleanspeech = [90; 0];
ang_laughter = [20; 0];

fs = 16000;
collector = phased.WidebandCollector('Sensor',ula,'PropagationSpeed',c,...
    'SampleRate',fs,'NumSubbands',1000,'ModulatedInput', false,'Wavefront','Plane');

t_duration = 3;  % 3 seconds
t = 0:1/fs:t_duration-1/fs;

prevS = rng(2008);
noisePwr = 1e-4; % noise power

% preallocate
NSampPerFrame = 1000;
NTSample = t_duration*fs;
sigArray = zeros(NTSample,Nele);

voice_cleanspeech = zeros(NTSample,1);
voice_laugh = zeros(NTSample,1);

% set up audio device writer
audioWriter = audioDeviceWriter('SampleRate',fs, ...
        'SupportVariableSizeInput', true);
%isAudioSupported = (length(getAudioDevices(audioWriter))>1);


speechFileReader = dsp.AudioFileReader('clean_speech.wav',...
    'SamplesPerFrame',NSampPerFrame);
laughterFileReader = dsp.AudioFileReader('笑声.wav',...
    'SamplesPerFrame',NSampPerFrame);

%% simulate
for m = 1:NSampPerFrame:NTSample
    sig_idx = m:m+NSampPerFrame-1;
    x2 = speechFileReader();
    x3 = 2*laughterFileReader();
    % 此处说明了共有三种声音：语音，笑声，噪声
    temp = collector([x2 x3],...
        [ ang_cleanspeech ang_laughter]) + ...
        sqrt(noisePwr)*randn(NSampPerFrame,Nele);
   % if isAudioSupported
   %     play(audioWriter,0.5*temp(:,3));
   % end
    sigArray(sig_idx,:) = temp;
   
    voice_cleanspeech(sig_idx) = x2;
    voice_laugh(sig_idx) = x3;
end

% audiowrite('mix.wav',sigArray(:,1),16000);


%==================以上模拟完毕，6路语音存放在sigArray中============================

%% 获取输入信号 做这步纯属是为了方便理解
x1=sigArray(:,1);
x2=sigArray(:,2);
x3=sigArray(:,3);
x4=sigArray(:,4);
x5=sigArray(:,5);
x6=sigArray(:,6);
X=[x1 x2 x3 x4 x5 x6].'; 
% audiowrite('mix.wav',x1,16000);

%% 初始化
c=340;
desire_angle=0; % 期望角度   这里有个细节，这里的0度对应模拟ULA的90度
NFFT=1024; % 帧长1024
Overlap=NFFT/2;
MicNum=6; % 6个麦克风
wind=hanning(NFFT); % 窗函数使用汗宁窗
MArray=0.05*[-2.5,0;-1.5,0;-0.5,0;0.5,0;1.5,0;2.5,0]; % 建立坐标系
W=zeros(MicNum,NFFT/2+1); % 初始化权重 频域
fs=16000; % 采样率
Freband=linspace(0,fs/2,NFFT/2+1);
As=zeros(MicNum,NFFT/2+1); % steervector  
n_blocks=fix((size(sigArray,1)-Overlap)/Overlap);
xbuf=zeros(MicNum,NFFT);
Xbuf=zeros(MicNum,NFFT/2+1); % fft后的xbuf
P=zeros(MicNum,MicNum);
F=zeros(MicNum,NFFT/2+1);
ffy=zeros(NFFT/2+1,1);
y_time=zeros(NFFT,1);
Output=zeros(Overlap*(n_blocks+1),1);
outbuf=zeros(Overlap,1);

mu=0.008;

%% 初始化As也就是给出约束条件
for freIdx=1:NFFT/2+1
    coefficient = 2 * pi * Freband(freIdx)/c;
    h=[cosd(desire_angle),sind(desire_angle)].';
    tao = MArray*h;
    As(:,freIdx)=exp(j*coefficient*tao);
end

P=eye(MicNum,MicNum)-As*As'./(norm(As,2)^2); %6*6
F=As./(norm(As,2)^2); % 6*513
W=F; % 6*513
num=1;
mydisplay=[];
for n=1:n_blocks
    xbuf(:,Overlap+1:end)=X(:,(n-1)*Overlap+1:n*Overlap);
    
    for m=1:MicNum
        tmp=fft(xbuf(m,:).');
        Xbuf(m,:)=tmp(1:NFFT/2+1).';
    end
    fprintf('%d\n',num);
    for freIdx=1:NFFT/2+1
        ffy(freIdx,1)=W(:,freIdx)'*Xbuf(:,freIdx)./MicNum; % 1*6 * 6*1 
        W(:,freIdx)=P*[W(:,freIdx)-mu*Xbuf(:,freIdx).*conj(ffy(freIdx,1))]+F(:,freIdx);      
    end
    
    y_time=real(ifft([ffy;flipud(conj(ffy(2:end-1)))]));
    
    Output((n-1)*Overlap+1:n*Overlap,1)=outbuf+y_time(1:Overlap);
    outbuf=y_time(Overlap+1:end);
    xbuf(:,1:Overlap)=xbuf(:,Overlap+1:end);
    num=num+1;
end


audiowrite('result0.wav',Output,fs);





