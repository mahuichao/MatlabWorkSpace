function [alpha] = forwardAlg(K,Y,a,b,phi)
%FORWARDALG 前向算法 已知模型的情况下，求解观测序列最大概率P(O|lamada)
%   a:转移概率
%   b:K*N
%   phi:初始概率分布
[N,~]=size(Y);
alpha=zeros(K,N); % 前向概率
threhold=10e-4;

for k=1:K
        alpha(k,1)=phi(k,1)*b(k,1);
end


for n=2:N
    for j=1:K
        tmp=0;
        for i=1:K
            tmp=tmp+alpha(i,n-1)*a(i,j);
        end
        alpha(j,n)=tmp*b(j,n);
        if alpha(j,n)<threhold
           alpha(j,n)=threhold;
        end
    end
end

end

