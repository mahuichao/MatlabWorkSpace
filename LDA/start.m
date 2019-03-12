close all; clc; clear all;

Data=load("iris.mat");
inputs=Data.Inputs;
targets=Data.Targets;
% 获取label
L=vec2ind(targets).';
[Y, W, lambda] = LDA(inputs,L);

figure;

D = size(inputs,1);
for d=1:D
    % 原始数据
    subplot(D,2,2*d-1);
    plot(inputs(d,:));
    ylabel(['x_' num2str(d)]);
    if d==D
        xlabel('Sample Index');
    end
    if d==1
        title('Original Data');
    end
    grid on;
    
    % 转换后数据
    subplot(D,2,2*d);
    plot(Y(d,:));
    ylabel(['y_' num2str(d)]);
    if d==D
        xlabel('Sample Index');
    end
    if d==1
        title('LDA Output');
    end
    grid on;
    
end
