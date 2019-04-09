function [output] = genRandMatx(m)
%GENRANDMATX 生成一个随机矩阵，不过每列相加为1
output=rand(m);
for i=1:m
    output(i,:)=output(i,:)./sum(output(i,:));
end

end

