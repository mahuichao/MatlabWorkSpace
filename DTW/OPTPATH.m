function [path] = OPTPATH(dtw)
%OPTPATH ×îÓÅÂ·¾¶
path=[];
[i,j]=size(dtw);
row=1;
while(i>1 && j>1)
    if(i==1)
        j=j-1;
    elseif(j==1)
        i=i-1;
    else
        if(dtw(i-1,j)==min(dtw(i-1,j),min(dtw(i-1,j-1),dtw(i,j-1))))
            i=i-1;
        elseif(dtw(i-1,j-1)==min(dtw(i-1,j),min(dtw(i-1,j-1),dtw(i,j-1))))
            i=i-1;j=j-1;
        else
            j=j-1;
        end
        path(row,:)=[i,j];
        row=row+1;
    end
    
   
end

end

