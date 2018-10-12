function [Output] = mfcc(data,framelength,keepedCoeffs,filterNums,lowFreq,highFreq,fs,time)
%MFCC 此处显示有关此函数的摘要
%   此处显示详细说明
% 我们每次只提取time的时间
data=data(1:fs*time,1);
% data=data./max(abs(data)); % 进行一次归一化
% data=data-mean(data); % 消除直流分量
% 帧长度
% framelength=512;
% 重叠区间长度
overlap=framelength/2;
% 最终留下的系数个数
%keepedCoeffs=12;
% 数据长度,总长度
datalength=size(data,1);
% 块个数，不够补0
Nblock=ceil((datalength-overlap)/overlap);
% DFT使用的窗，汉明窗
win=hamming(framelength);
% 存放帧
xbuf=zeros(framelength,1);
% DFT后帧长度，由于傅立叶变换后会有复数共轭的特性，所以只取一半+1个，具体查看傅立叶变换
Xbuf=zeros(framelength/2+1,1);
% 功率谱
P=zeros(framelength/2+1,1);
% 下截止频率
%lowFreq=300;
% 上截止频率
%highFreq=fs/2;
% 滤波器组内滤波器个数
%filterNums=26; 
% 需要filterNums个频带，那么至少在两边加两个点 也就是28个点
points2Keeped=filterNums+2;
% 中心频率集合
fbins=zeros(points2Keeped,1);
% 取log之后的梅尔功率谱
melPower=zeros(filterNums,1);
% 最终的输出
Output=zeros(Nblock,keepedCoeffs);
for i =1:Nblock
    % 分帧操作
    xbuf(framelength/2+1:end)=data(framelength/2*(i-1)+1:i*framelength/2,1);
    % 傅立叶变换+窗
    tmp=fft(xbuf(:,1).'.*win.');
    % 把前半+1存入Xbuf中
    Xbuf(:,1)=tmp(1:framelength/2+1);
    % 功率谱
    P(:,1)=power(abs(Xbuf(:,1)),2)./(framelength/2+1);
    % 梅尔映射
    mels=mfcc_melmap((lowFreq:highFreq).');
    m=linspace(mels(1),mels(end),points2Keeped);
    % 梅尔反映射回到标准频率
    h=mfcc_inv_melmap(m.');
    % 确定与我们要求的频率最接近的bin
    fbins=floor(framelength/fs.*h+(h./fs));
    % 滤波器组
    fbanks=mfcc_filterBanks(filterNums,fbins,framelength);%.*(fs/framelength);
%     for x=1:26
%         plot((1:framelength/2+1).*(fs/framelength),abs(fbanks(x,:)));
%         hold on;
%     end
    % 对功率谱应用滤波器组
    for n=1:filterNums
        melPower(n,1)=log10(fbanks(n,:)*P(:,1));
    end
    % 离散余弦变换
    mels=mfcc_dct2(melPower,filterNums,keepedCoeffs);
    % 把结果存入最终输入中
    Output(i,:)=mels.';
    xbuf(1:framelength/2)=xbuf(framelength/2+1:end);
end

% Output=abs(Output);

end

