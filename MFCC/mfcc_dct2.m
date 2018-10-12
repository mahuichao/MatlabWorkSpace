function [output] = mfcc_dct2(x,N,keepedNums)
   output=zeros(keepedNums,1);
   tmpX=zeros(N,1);
   for i=1:N
       for j=1:keepedNums
          tmpX(i,1)=tmpX(i,1)+x(j,1)*cos(pi/N*(j-1/2)*(i-1));
       end
   end
   
   output(:,1)=tmpX(2:keepedNums+1,1);
end

