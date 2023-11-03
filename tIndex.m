function [value] = tIndex(metric,ind1,ind2,ind3)
    mode = 0;
    
    % Hard coded to 4. Maybe improve in the future?
    ind1vec = zeros(1,4);
    ind2vec = zeros(4,1);
    ind3vec = zeros(4,1);

    if ind1 == ':'
        mode = 1;
    else
        ind1vec(ind1) = 1;
    end

    if ind2 == ':'
        mode = 2;
    else
        ind2vec(ind2) = 1;
    end

    if nargin == 4
        if ind3 == ':'
            first = ind1vec * tensorprod(metric,[1; 0; 0],3,1) * ind2vec;
            second = ind1vec * tensorprod(metric,[0; 1; 0],3,1) * ind2vec;
            third = ind1vec * tensorprod(metric,[0; 0; 1],3,1) * ind2vec;
            value = [first second third];
        else
            ind3vec(ind3) = 1;
            metric = tensorprod(metric,ind3vec,3,1);
            value = ind1vec * metric * ind2vec;
        end
    else
        if mode == 0
            value = ind1vec * metric * ind2vec;
        elseif mode == 1
            value = metric * ind2vec;
        elseif mode == 2
            value = ind1vec * metric;
        end
    end
end