function [m,idx]=rollmax(w,k)
    n = length(w);
    if (k>=n)||(k<=1)
        [m,idx]=max(w);
    else
        roll=w(1:end-k+1);
        for i=2:k
            roll=roll+w(i:end-k+i);
        end
        [~,idx]=max(roll);
        idx=idx-1;
        if idx<1
            idx=1;
        end
        m=w(idx);
    end
end