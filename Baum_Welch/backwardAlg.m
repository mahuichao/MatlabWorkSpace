function [bita] = backwardAlg(K,Y,a,b)
%BACKWARDALG 后向算法
%   a:转移概率
%   b:K*N
%   phi:初始概率分布
[N,~]=size(Y);
bita=zeros(K,N);
bita(:,N)=1; % 后向概率
threhold=10e-4;

for n=N-1:-1:1
    for i=1:K
        tmp=0;
        for j=1:K
            tmp=tmp+a(i,j)*b(j,n+1)*bita(j,n+1);
        end
        bita(j,n)=tmp;
        if bita(j,n)<threhold
           bita(j,n)=threhold;
        end
    end
end

end

