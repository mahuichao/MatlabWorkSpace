function [output] = mfcc_melmap(frequencys)
    %  ‰»ÎŒ™––
    output=zeros(size(frequencys,1),1);
    for i=1:size(frequencys,1)
        output(i,1)=1125*log(frequencys(i)/700+1);
    end
end

