function [output] = mfcc_inv_melmap(melFreqs)
datalength=size(melFreqs,1);
output=zeros(datalength,1);
for i=1:datalength
    output(i,1)=700*(exp(melFreqs(i,1)/1125)-1);
end
end

