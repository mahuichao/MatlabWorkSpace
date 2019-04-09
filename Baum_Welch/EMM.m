function [alpha,mu,sigma] = EMM(Y,M,T)
%EMM EM+GMM 同时GMM为多变量混合高斯模型
%为了简化程序，初版我们不使用kmeans方法进行初始化
%参数说明：
% Y:数据 N*D
% M:模式个数，在分类问题中，表明类别个数
% T:迭代次数

%%初始化过程
[N,D]=size(Y);
alpha=ones(M,1).*1/M; % 各个模型的初始概率
mu=rand(M,D); % 初始均值，列为维度，行为不同的模型 M*D
sigma=zeros(D,D,M); 
for m=1:M
    sigma(:,:,m)=eye(D);
end

%% EM阶段
for t=1:T % 迭代
    p_m_all=zeros(N,1); % k观测数据属于各模型的总概率
    mu_up=zeros(M,D); % 均值更新公式上部分
    sigma_up=zeros(D,D,M); % 协方差更新公式上部分
    comm=zeros(M,D); % 需要重复计算的部分单独提出来
    alpha_up=zeros(M,1); % alpha更新公式的上部分
   
    for n=1:N
        for m=1:M      
            p_m_all(n,1)=p_m_all(n,1)+alpha(m,1)*mvnpdf(Y(n,:),mu(m,:),sigma(:,:,m));
        end
    end
    %% 迭代alpha
    for m=1:M
        for n=1:N
            comm(m,n)=(alpha(m,1)*mvnpdf(Y(n,:),mu(m,:),sigma(:,:,m)))./p_m_all(n,1);
            alpha_up(m,1)=alpha_up(m,1)+comm(m,n);
            mu_up(m,:)=mu_up(m,:)+Y(n,:).*comm(m,n);
        end    
      
        mu(m,:)=mu_up(m,:)./alpha_up(m,1);
  
        for n=1:N
             sigma_up(:,:,m)=sigma_up(:,:,m)+(Y(n,:)-mu(m,:)).'*(Y(n,:)-mu(m,:))*comm(m,n);
        end
        sigma(:,:,m)=sigma_up(:,:,m)./alpha_up(m,1);
       
        alpha(m,1)=alpha_up(m,1)/N;
    end
end

end

