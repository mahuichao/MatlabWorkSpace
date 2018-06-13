close all; clc; clear all;

%% 6个麦克风，5cm。
microphone = ...
    phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 20e3]);

Nele = 6;
ula = phased.ULA(Nele,0.05,'Element',microphone);
c = 340;                 

% ang_dft = [90; 0];
ang_cleanspeech = [30; 0];
% ang_laughter = [-20; 0];


fs = 16000;
collector = phased.WidebandCollector('Sensor',ula,'PropagationSpeed',c,...
    'SampleRate',fs,'NumSubbands',1000,'ModulatedInput', false); % 'Wavefront','Plane'

t_duration = 3;  % 3 seconds
t = 0:1/fs:t_duration-1/fs;

prevS = rng(2008);
noisePwr = 1e-4; % noise power


% preallocate
NSampPerFrame = 1000;
NTSample = t_duration*fs;
sigArray = zeros(NTSample,Nele);
% voice_dft = zeros(NTSample,1);
voice_cleanspeech = zeros(NTSample,1);
% voice_laugh = zeros(NTSample,1);

% set up audio device writer
% audioWriter = audioDeviceWriter('SampleRate',fs, ...
%         'SupportVariableSizeInput', true);
%isAudioSupported = (length(getAudioDevices(audioWriter))>1);

% dftFileReader = dsp.AudioFileReader('dft_voice_8kHz.wav',...
%     'SamplesPerFrame',NSampPerFrame);
% speechFileReader = dsp.AudioFileReader('cleanspeech_voice_8kHz.wav',...
%     'SamplesPerFrame',NSampPerFrame);
% laughterFileReader = dsp.AudioFileReader('laughter_8kHz.wav',...
%     'SamplesPerFrame',NSampPerFrame);
% dftFileReader = dsp.AudioFileReader('PCM005_R.wav',...
%     'SamplesPerFrame',NSampPerFrame);
speechFileReader = dsp.AudioFileReader('clean_speech_reb.wav',...
    'SamplesPerFrame',NSampPerFrame);
% laughterFileReader = dsp.AudioFileReader('laugher_new.wav',...
%     'SamplesPerFrame',NSampPerFrame);



%% simulate
for m = 1:NSampPerFrame:NTSample
    sig_idx = m:m+NSampPerFrame-1;
%     x1 = dftFileReader();
    x2 = speechFileReader();
%     x3 = laughterFileReader();
%     x3=x3(:,1);
%     x2=x2(:,1);
% temp = collector([x1 x2 x3],...
%         [ang_dft ang_cleanspeech ang_laughter])+sqrt(noisePwr)*randn(NSampPerFrame,Nele); 
% temp = collector([x2 x3],...
%         [ang_cleanspeech ang_laughter]);%+sqrt(noisePwr)*randn(NSampPerFrame,Nele); 
  temp = collector([x2 ],...
        [ang_cleanspeech ]);%+sqrt(noisePwr)*randn(NSampPerFrame,Nele); 
     sigArray(sig_idx,:) = temp;
     
%      voice_dft(sig_idx) = x1;
     voice_cleanspeech(sig_idx) = x2;
%      voice_laugh(sig_idx) = x3;
end

% audiowrite('mic_01.wav',sigArray(:,1),fs);
% audiowrite('mic_02.wav',sigArray(:,2),fs);
% audiowrite('mic_03.wav',sigArray(:,3),fs);
% audiowrite('mic_04.wav',sigArray(:,4),fs);
% audiowrite('mic_05.wav',sigArray(:,5),fs);
% audiowrite('mic_06.wav',sigArray(:,6),fs);

 
x1=sigArray(:,1);
x2=sigArray(:,2);
x3=sigArray(:,3);
x4=sigArray(:,4);
x5=sigArray(:,5);
x6=sigArray(:,6);

NFFT=512;
Overlap=256;
angle=120;
forgottor=0.4;
mu=0.05;
X=[x1 x2 x3 x4 x5 x6]';



%% 自己的算法 TF-GSC

[FBFOutput,Output]=TF_GSC(X,fs,angle,NFFT,mu,forgottor);
FBFOutput=FBFOutput(257:length(voice_cleanspeech));
Output=Output(257:length(voice_cleanspeech)); 
% voice_cleanspeech=voice_cleanspeech(1:end-256); 
% voice_laugh=voice_laugh(1:end-256);
% voice_dft=voice_dft(1:length(Output));
% fbf_Cbf=calculate_SNR(sigArray,voice_cleanspeech,noisePwr,FBFOutput,voice_laugh); 
% agCbf=calculate_SNR(sigArray,voice_cleanspeech,noisePwr,Output,voice_laugh); 
audiowrite('Hinput_33_withTD_inH2.wav',sigArray(:,1),16000);
audiowrite('Houtput_33_withTD_inH2.wav',Output,16000);


%% 使用内置TimeDelay进行滤波
angSteer=[20;0];%ang_cleanspeech;
beamformer=phased.TimeDelayBeamformer('SensorArray',ula,...
    'SampleRate',fs,'Direction',angSteer,'PropagationSpeed',c);
signalsource=dsp.SignalSource('Signal',sigArray,...
    'SamplesPerFrame',NSampPerFrame);
cbfOut=zeros(NTSample,1);
for m=1:NSampPerFrame:NTSample
    temp=beamformer(signalsource());
    cbfOut(m:m+NSampPerFrame-1,:)=temp;
end

% agCbf=calculate_SNR(sigArray,voice_cleanspeech,noisePwr,cbfOut,voice_laugh);
% audiowrite('TD.wav',cbfOut,16000);


%% frost算法
frostbeamformer = ...
    phased.FrostBeamformer('SensorArray',ula,'SampleRate',16000,...
    'PropagationSpeed',c,'FilterLength',5,'DirectionSource','Input port');

reset(signalsource);
signalsource=dsp.SignalSource('Signal',sigArray,...
    'SamplesPerFrame',NSampPerFrame);
FrostOut = zeros(NTSample,1);
for m = 1:NSampPerFrame:NTSample
    FrostOut(m:m+NSampPerFrame-1,:) = ...
        frostbeamformer(signalsource(),ang_cleanspeech);
end

% agCbf=calculate_SNR(sigArray,voice_cleanspeech,noisePwr,FrostOut,voice_laugh);

release(frostbeamformer);
ang_cleanspeech_est = [-5; 5];  % Estimated steering direction

reset(signalsource);
FrostOut2 = zeros(NTSample,1);
for m = 1:NSampPerFrame:NTSample
    FrostOut2(m:m+NSampPerFrame-1,:) = frostbeamformer(signalsource(),...
        ang_cleanspeech_est);
end

% calculate_SNR(sigArray,voice_cleanspeech,noisePwr,FrostOut,voice_laugh);

% audiowrite('mic_01.wav',sigArray(:,1),fs);
% audiowrite('mic_02.wav',sigArray(:,2),fs);
% audiowrite('mic_03.wav',sigArray(:,3),fs);
% audiowrite('mic_04.wav',sigArray(:,4),fs);
% audiowrite('mic_05.wav',sigArray(:,5),fs);
% audiowrite('mic_06.wav',sigArray(:,6),fs);