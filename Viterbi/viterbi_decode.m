function [delta,kx] = viterbi_decode(a,b,phi)
%VITERBI_DECODE 维特比算法解码
%   a:K*K 状态转移概率矩阵
%   b:K*N 发射概率矩阵
%   phi:K*1 初始概率分布
[K,N]=size(b);

delta=zeros(N,K);
kx=zeros(N,1);

for k=1:K
    delta(1,k)=phi(k,1)*b(k,1);
    kx(1,1)=1;
end

for n=2:N    
   [delta(n,:),kx(n,1)]=max(delta(n-1,:)*a,[],2);
   delta(n,:)=delta(n,:).'.*b(:,n);
end

end

