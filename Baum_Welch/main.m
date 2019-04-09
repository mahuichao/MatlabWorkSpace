close all; clc; clear all; 

%% 生成数据
R1=mvnrnd([1,2.25],[1.75,1;1,1.75],30);
R2=mvnrnd([12.9,10.3],[3.2,1.25;1.25,3.2],70);

Y=[R1;R2];
plot(R1(:,1),R1(:,2),'ro'); hold on;
plot(R2(:,1),R2(:,2),'bo');hold on;

%% 暂时假设一个状态对应一个高斯分布
K=2; % 多高斯模型，此处为双高斯模型
T=100; % 迭代次数


% 我们认为得到的mu和sigma就是HMM的Baum-Welch算法中的发射矩阵b
% 那么他们的对应关系是什么呢？
[a,b,alpha]=BaumWelch(K,Y,T);
