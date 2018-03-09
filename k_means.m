function [ newX ] = k_means( x,c )
%K_MEANS k-means均值聚类 我们这里使用马氏距离进行计算
%   x:数据 多维
%   c:类别的个数;
%   iter:迭代次数(考虑到中心点一直不收敛,暂时不使用)
[charcs,datalength]=size(x);
U=zeros(c,charcs); % c类，每类charcs个特征
LastU=zeros(c,charcs);
randIdx=randperm(datalength);
drawIdx=randIdx(1:c);
newX=[x;zeros(1,datalength)]; % 第三行表示类别
for cIdx=1:c
    U(cIdx,:)=x(:,drawIdx(cIdx)).';
end
flag=0;
LastU=U;
canChange=1;
while flag~=1
    
    % 循环分类
    for n=1:datalength
        dOld=Inf;
        for cIdx=1:c
            dNew=norm(x(:,n).'-U(cIdx,:)).^2;
            if dOld>dNew % 离新的较近
                newX(:,n)=[x(:,n);cIdx];
                dOld=dNew;
            else
                continue;
            end
        end
    end
    
%     newXWid=size(newX,1);
%     newXLen=size(newX,2);
    % 对每个类进行分析
    for cIdx=1:c
        % 拿出各个源的数据
        classIdxs=find(newX(3,:)==cIdx);
        dataSet=zeros(charcs,length(classIdxs));
        dataSet(:,:)=newX(1:charcs,classIdxs);
        
        % 更新重心
        U_up=zeros(charcs,1);
        U_down=0;
        for len=1:length(dataSet)
            U_up=U_up+dataSet(:,len);
            U_down=U_down+1;
        end
        U(cIdx,:)=U_up.'./U_down;
        if LastU(cIdx,:)-U(cIdx,:)==zeros(1,charcs)
            if canChange==1
                flag=1;
            end
        else
            flag=0;
            canChange=0;
        end
    end
    LastU=U;
    canChange=1;
end


end

