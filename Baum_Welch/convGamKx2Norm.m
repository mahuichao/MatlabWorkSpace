function [gamma,kx] = convGamKx2Norm(logGamma,logKx)
%CONVGAMKX2NORM 把logGamma和logKx去掉log
%   logGamma:K*N
%   logKx:K*K*N

[K,N] = size(logGamma);
gamma=zeros(K,N); % 给定模型和观测序列，在时刻t处于某个状态的概率 普通版
kx=zeros(K,K,N); % 给定模型和观测序列，在时刻t处于某个状态，在时刻t+1时刻处于某个时刻的概率 普通版

for k=1:K
    maxGamma=0;
    maxKx=0;
    for n=1:N
        maxGamma=max(logGamma(k,:));
        maxKx=max(reshape(logKx(k,:,:), K*(N-1), 1));
    end
    logGamma(k,:)=logGamma(k,:)-maxGamma;
    logKx(k,:)=logKx(k,:)-maxKx;
end

gamma=exp(logGamma);
kx=exp(logKx);

end

