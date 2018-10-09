function [mu,sigma] = EMM(Y,K,times)
%EMM 此处显示有关此函数的摘要
%   此处显示详细说明
% y:训练数据
% k:模型个数
% times:迭代次数

%% 参数初始化
% 主要有三个：mu，sigma，alpha
[M,N]=size(Y); % M 表示数据个数，N 表示维数
alpha=ones(K,1).*1/K; 
mu=rand(K,N);
sigma=eye(N).*K; 

resp=zeros(M,K);

%% 进入迭代
for t=1:times
    %% E步骤，求取响应度
    for k=1:K
        resp(:,k)=alpha(k,1).*getPDF(Y,mu(k,:),sigma);
    end
    for i=1:N
        resp(i,:)=resp(i,:)./sum(resp(i,:));
    end
    
    %% M步骤，计算新一轮迭代的模型参数
    for k=1:K
      
       r_sum=sum(resp(:,k));
       mu(k,:)= sum(resp(:,k).'*Y(:,:))./r_sum;
       
       cov_k=zeros(N,N);
       for i=1:M
            cov_k=cov_k+resp(i,k)*((Y(i,:)-mu(k,:))'*(Y(i,:)-mu(k,:)))./r_sum;
       end
       sigma=cov_k;
       alpha(k,1)=r_sum/M;
    
    end
end

end

