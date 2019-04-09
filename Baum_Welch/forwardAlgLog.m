function [alpha,logc] = forwardAlgLog(K,Y,a,b,phi)
%FORWARDALG 前向算法 已知模型的情况下，求解观测序列最大概率P(O|lamada)
% b:K*N
% a:K*K
% phi:K*1
[N,~]=size(Y);
alpha=zeros(K,N); % 前向概率
logc = zeros(N,1);


tmp=log(phi(1,:))+log(b(1,:));
logc(1) = log( sum( exp( tmp - max(tmp) ) ) ) + max(tmp);
alpha(:,1)= -logc(1) + log(b(:,1)) + log(phi);


for n=2:N
    tmp = bsxfun(@plus, bsxfun(@plus, log(a), alpha(:,n-1)'), log(b(:,n)));
    logc(n) = log ( sum( sum ( exp ( tmp - max(tmp(:)) ) ) ) ) + max(tmp(:));
    for j=1:K
        tmp2 = alpha(:,n-1) + log(a(:,j));
        if (isinf(max(tmp2)))
                alpha(j,n) = -inf;
        else
                alpha(j,n) = -logc(n) + log(b(j,n)) + log( sum( exp( tmp2 - max(tmp2) ) ) ) + max(tmp2);
        end
    end
end


end

