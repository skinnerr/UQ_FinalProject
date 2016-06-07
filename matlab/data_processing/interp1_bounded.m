function [ vq ] = interp1_bounded( x, v, xq )
% Same as MATLAB's interp1 function, but restricts xq to fall within range of x.

    maxx = max(x);
    minx = min(x);
    ind = find(xq > maxx);
    if ind
        xq(ind) = maxx;
    end
    ind = find(xq < minx);
    if ind
        xq(ind) = minx;
    end
    
    vq = interp1(x, v, xq);

end

