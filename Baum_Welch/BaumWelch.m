function [a,b,alpha] = BaumWelch(K,Y,T)
%BAUMWELCH HMM训练大法之Baum-Welch
%   mu:K*D,sigma:D*D*K,alpha:K*1
%   K:状态的个数
%   mu:高斯混合模型的均值
%   sigma:高斯混合模型的协方差矩阵
%   alpha:高斯混合模型的初始概率分布
%   times:迭代次数

[N,D]=size(Y);
times=200; % GMM迭代次数
% GMM结果并无问题
%  [~,mu,sigma] = Gmm(Y, K, 'diag');
[alpha,mu,sigma]=EMM(Y,K,times); % mu:K*D,sigma:D*D*K,alpha:K*1
for k=1:K
    sigma(:,:,k)=diag(diag(sigma(:,:,k)));
end

 
% color={'r','b'};
% for m=1:K
%     error_ellipse(sigma(:,:,m), mu(m,:)', 'style', color{m}); hold on 
% end

%% 初始化转移矩阵a，发射矩阵b，初始概率分布alpha
a=genRandMatx(K);
b=zeros(K,N); %现在的模型是一个类对应一个单高斯
for k=1:K
    for n=1:N
        b(k,n)=mvnpdf(Y(n,:),mu(k,:),sigma(:,:,k));
    end
end

% c=zeros(K,N);
% b = Gauss_logp_xn_given_zn(Y,mu.',sigma);

%% Baum-Welch算法阶段
for t=1:T
    
    % 1.前向算法 初始计算出来的pObs数值存在问题，等后面开始迭代后再观察
    [forPro]=forwardAlg(K,Y,a,b,alpha); % forPro:K*N pObs_1:1*1
%     [forPro,logc]=forwardAlgLog(K,Y,a,b,alpha);

    % 2.后向算法
    [bakPro]=backwardAlg(K,Y,a,b); % forPro:K*N p_obs_2:1*1 
%     [bakPro]=backwardAlgLog(K,Y,a,b,logc);

    % 3.计算这里暂且称之为过渡概率
%     logGamma=zeros(K,N); % 给定模型和观测序列，在时刻t处于某个状态的概率 log版
%     logKx=zeros(K,K,N); % 给定模型和观测序列，在时刻t处于某个状态，在时刻t+1时刻处于某个时刻的概率 log版
    gamma=zeros(K,N);
    kx=zeros(K,N);
    
    % calculate loggamma
%     logGamma = forPro + bakPro;
    
    % calculate logksi
%     for n = 2:N
%         logKx(:,:,n) = -logc(n) + bsxfun(@plus, bsxfun(@plus, log(a), forPro(:,n-1)), log(b(:,n)) + bakPro(:,n));
%     end
%     logKx(:,:,1) = [];
%     
%     [gamma,kx]=convGamKx2Norm(logGamma,logKx);

    com=zeros(N,1); % 通用计算部分
    % 这里是用最基础的乘法完成的，后期可以改成矩阵形式
    for n=1:N-1
        for i=1:K
            for j=1:K
                com(n,1)=com(n,1)+forPro(i,n)*a(i,j)*mvnpdf(Y(n+1,:),mu(j,:),sigma(:,:,j))*bakPro(j,n+1);
            end
        end
    end
    
    % 求克西 
    for i=1:K
        for j=1:K
            for n=1:N-1
                kx(i,j,n)=forPro(i,n)*a(i,j)*mvnpdf(Y(n+1,:),mu(j,:),sigma(:,:,j))*bakPro(j,n+1)./com(n,1);
            end
        end
    end
    
    % 求伽马
    for i=1:K
        for n=1:N-1
            gamma(i,n)=sum(kx(i,:,n));
        end
    end
    gamma(K,N)=forPro(K,N)/sum(forPro(:,N).');
    
%     % 开始迭代操作
%     % 1) 初始概率分布迭代
    for i=1:K
        alpha(i,1)=gamma(i,1);
    end
%     
%     
%     % 2) 状态转移概率迭代
    for i=1:K
        for j=1:K
            a(i,j)=sum(kx(i,j,:))./sum(gamma(i,:));
        end
    end
%     
% 
%     % 3) 发射概率矩阵迭代
    for j=1:K
        for n=1:N
            b(j,n)=gamma(j,n)/sum(gamma(j,:)); % 所以这里的更新存在问题。
        end
    end
    
end

end

