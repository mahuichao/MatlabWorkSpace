close all; clc; clear all;

%% 生成数据 目前数据是两个维度
R1=mvnrnd([1,2.25],[1.75,1;1,1.75],30);
R2=mvnrnd([12.9,10.3],[3.2,1.25;1.25,3.2],70);

Y=[R1;R2];
plot(R1(:,1),R1(:,2),'ro'); hold on;
plot(R2(:,1),R2(:,2),'bo');hold on;
% plot(Y(:,1),Y(:,2),'r+');hold on;

M=2; % 多高斯模型，此处为双高斯模型
T=300; % 迭代次数
[alpha,mu,sigma]=EMM(Y,M,T);
color={'r','b'};
for m=1:M
    error_ellipse(sigma(:,:,m), mu(m,:)', 'style', color{m}); hold on
end