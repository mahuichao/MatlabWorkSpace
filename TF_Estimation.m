function [ H ] = TF_Estimation(x,NFFT,segments,window)
% system identification
% x: 输入信号 6*13000 此函数内部不做缓冲
[MicNum,DataLength]=size(x);

m1=zeros(NFFT,segments);
m2=zeros(NFFT,segments);
m3=zeros(NFFT,segments);
m4=zeros(NFFT,segments);
m5=zeros(NFFT,segments);
m6=zeros(NFFT,segments);
H=zeros(MicNum,NFFT);


for i=1:segments
    m1(:,i)=x(1,(i-1)*NFFT+1:i*NFFT);
    m2(:,i)=x(2,(i-1)*NFFT+1:i*NFFT);
    m3(:,i)=x(3,(i-1)*NFFT+1:i*NFFT);
    m4(:,i)=x(4,(i-1)*NFFT+1:i*NFFT);
    m5(:,i)=x(5,(i-1)*NFFT+1:i*NFFT);
    m6(:,i)=x(6,(i-1)*NFFT+1:i*NFFT);
end


ftm1=zeros(NFFT,segments);
ftm2=zeros(NFFT,segments);
ftm3=zeros(NFFT,segments);
ftm4=zeros(NFFT,segments);
ftm5=zeros(NFFT,segments);
ftm6=zeros(NFFT,segments);

for i=1:segments
    ftm1(:,i)=fft(m1(:,i).*window);
    ftm2(:,i)=fft(m2(:,i).*window);
    ftm3(:,i)=fft(m3(:,i).*window);
    ftm4(:,i)=fft(m4(:,i).*window);
    ftm5(:,i)=fft(m5(:,i).*window);
    ftm6(:,i)=fft(m6(:,i).*window);
end

Rx11=mean(conj(ftm1).*ftm1,2);
Rx12=mean(conj(ftm1).*ftm2,2);
Rx13=mean(conj(ftm1).*ftm3,2);
Rx14=mean(conj(ftm1).*ftm4,2);
Rx15=mean(conj(ftm1).*ftm5,2);
Rx16=mean(conj(ftm1).*ftm6,2);

denotor=mean(conj(ftm1).*ftm1.*conj(ftm1).*ftm1,2)-mean(conj(ftm1).*ftm1,2).^2;

H(1,:)=1;
H(2,:)=(mean(conj(ftm1).*ftm1.*conj(ftm2).*ftm1,2)-Rx11.*Rx12)./denotor;
H(3,:)=(mean(conj(ftm1).*ftm1.*conj(ftm3).*ftm1,2)-Rx11.*Rx13)./denotor;
H(4,:)=(mean(conj(ftm1).*ftm1.*conj(ftm4).*ftm1,2)-Rx11.*Rx14)./denotor;
H(5,:)=(mean(conj(ftm1).*ftm1.*conj(ftm5).*ftm1,2)-Rx11.*Rx15)./denotor;
H(6,:)=(mean(conj(ftm1).*ftm1.*conj(ftm6).*ftm1,2)-Rx11.*Rx16)./denotor;



