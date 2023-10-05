function [g2,gtime,g2zero,g2zero_corr,g2_err]=g2_pulseresolved(tPh0,tPh1,Sync0,Sync1,tauRes,RepTime,opts)
    arguments
        tPh0;
        tPh1;
        Sync0;
        Sync1;
        tauRes;
        RepTime;
        opts.NumPks {mustBeInteger,mustBeGreaterThan(opts.NumPks,3),mustBeLessThan(opts.NumPks,15)} = 5;
        opts.pBG {mustBeNonnegative} = 0;
    end
    n=fix((opts.NumPks-1)/2);
    i0=length(Sync0);
    i1=length(Sync1);
    nbins=fix(2*RepTime/tauRes);
    gtime=(-(n+1)*RepTime:tauRes:(n+1)*RepTime);
    g2=zeros(1,length(gtime));
    g2int=zeros(opts.NumPks,1);
    g2time=zeros(1,min(i0,i1));    
    Sync1=Sync1-n-1;
    for i=1:opts.NumPks
        Sync1=Sync1+1;
        if Sync0(1)<Sync1(1)
            pt=0;
        else
            pt=1;
        end
        j=0;    k0=1;    k1=1;
        while (k0<=i0)&&(k1<=i1)
            if pt==0
                if Sync0(k0)==Sync1(k1)
                    j=j+1;
                    g2time(j)=tPh1(k1)-tPh0(k0);
                    pt=1;
                    k1=k1+1;
                    k0=k0+1;
                elseif Sync0(k0)>Sync1(k1)
                    pt=1;
                    k1=k1+1;
                else
                    k0=k0+1;
                end
            else
                if Sync1(k1)==Sync0(k0)
                    j=j+1;
                    g2time(j)=tPh1(k1)-tPh0(k0);
                    pt=0;
                    k0=k0+1;
                    k1=k1+1;
                elseif Sync1(k1)>Sync0(k0)
                    pt=0;
                    k0=k0+1;
                else
                    k1=k1+1;
                end
            end
        end
        g2time=g2time(1:j);
        g2hist=histcounts(g2time,nbins,'BinLimits',[-RepTime,RepTime]);
        g2int(i)=trapz(g2hist);
        g2((i-1)*nbins/2+1:(i+1)*nbins/2)=g2((i-1)*nbins/2+1:(i+1)*nbins/2)+g2hist;
    end
    g2one=(sum(g2int)-g2int(n+1))/(n*2);
    if g2one==0
        g2_err=1;
        g2zero=1;
    else
        g2_err=round(1/g2one,3);
        g2zero=round(g2int(n+1)/g2one,3);
    end
    g2zero_corr=round(g2zero-4*opts.pBG,3);
    if g2zero_corr<0
        g2zero_corr=0;
    end
    fprintf(1,"g2(0) counts: %g, <g2(1) counts>: %g\ng2(0) = %g, g2(0)_corr = %g, g2_err = %g\n",g2int(n+1),g2one,g2zero,g2zero_corr,g2_err);
    
end