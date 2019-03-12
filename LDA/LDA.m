function [Y,W,lamada] = LDA(inputs,labels)
%LDA算法
%   参考实现：Yarpiz (www.yarpiz.com)

wholeData=[inputs;labels.'];
dlen=size(wholeData,2);
d=size(inputs,1);
classes=unique(labels);
classNumber=size(classes,1);
SB=zeros(d,d);
SW=zeros(d,d);

m=sum(wholeData(1:4,:),2)./dlen;
% 求SB和SW,求广义特征值,与博客公式相一致
for i=1:classNumber
    [~,n]=find(wholeData(end,:)==i);
    Xi=wholeData(1:4,n);
    ni=size(Xi,2);
    mi=sum(Xi,2)./ni;
    SB=SB+ni*(mi-m)*(mi-m)';
    for n=1:ni
        SW=SW+(Xi(:,n)-mi)*(Xi(:,n)-mi).';
    end
end
% 根据SB和SW进行广义特征值分解
[W,LAMBDA]=eig(SB,SW);
lamada=diag(LAMBDA);
[lamada,desorder]=sort(lamada,'descend');
W=W(:,desorder);
Y=W.'*inputs;
end

