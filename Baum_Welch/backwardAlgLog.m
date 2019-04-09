function [bita] = backwardAlgLog(K,Y,a,b,logc)
%BACKWARDALG 后向算法
%   a:转移概率 K*K
%   b:K*N
%   phi:初始概率分布
[N,~]=size(Y);
bita=zeros(K,N);
bita(:,N)=0; % 后向概率 因为log1=0


for n=N-1:-1:1
    for i=1:K
        tmp = bita(:,n+1) + log(b(:,n+1)) + log(a(:,i));
        bita(i,n) = -logc(n+1) + log( sum( exp( tmp - max(tmp) ) ) ) + max(tmp);
    end
end

end

