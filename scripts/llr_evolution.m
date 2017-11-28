function llr_evolution()
clear all;
close all;

log_plot = false;
dv=3;
dc=6;
M = 5;
I = 100;
ii_plot=4;
tol = 1e-6;
Pe_min = 1e-16;

Nq = 4; %
sig = .84;

p0 = get_quant_biawgn_pmf(sig.^2, Nq);
L0 = log(p0)- fliplr(log(p0));

p=p0;
L = L0;


if(1)
    Nq = 0;
    subplot(3,1,1);
    stem(L,p);
    if(log_plot)
        set(gca,'yscal','log');
    end

    % CN update
    p = get_product_pmf(dc,p, Nq);
    L = log(p)-log(fliplr(p));
    subplot(3,1,2);
    stem(L,p);
    if(log_plot)
        set(gca,'yscal','log');
    end


    % VN update
    p = get_product_pmf(dv,p, p0 , Nq);
    L = log(p)-log(fliplr(p));
    subplot(3,1,3);
    stem(L,p);
    if(log_plot)
        set(gca,'yscal','log');
    end
else
    for ii=1:I
        % CN update
        p = get_product_pmf(dc,p , Nq);
        p = p / sum(p);
        L = log(p)-log(fliplr(p));


        % VN update
        p = get_product_pmf(dv,p, p0 , Nq);
        p = p / sum(p);
        L = log(p)-log(fliplr(p));
        
        ii
        Pe = sum(p(1:(length(p)/2)))
        sum(p)
        
        
        
        %if(mod(ii-1,ii_plot)==0)
            figure(ii);
            stem(L,p);
            if(log_plot)
                set(gca,'yscal','log');
            end
        %end
        
        if (Pe <= Pe_min)
            break;
        end
    end
end

end
