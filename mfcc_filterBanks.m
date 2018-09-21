function [output] = mfcc_filterBanks(filterNums,f,NFFT)
    output=zeros(filterNums,NFFT/2+1);
    for m=2:filterNums+1
        tmp_left=floor(f(m-1));
        tmp_mid=floor(f(m));
        tmp_right=floor(f(m+1));
        
        for k=tmp_left:tmp_mid
            output(m-1,k)=(k-f(m-1))/(f(m)-f(m-1));
        end
        
        for k=tmp_mid:tmp_right
            output(m-1,k)=(f(m+1)-k)/(f(m+1)-f(m));
        end
    end

end

