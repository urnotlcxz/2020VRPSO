function outtime = formatTime(t)

if t < 60
    s = t;
    outtime = sprintf('%.2f s', s);
elseif t >= 60 && t < 3600
    s = mod(t, 60);
    m = floor((t-s)/60);
    outtime = sprintf('%d min %.2f s', m, s);
elseif t >= 3600
    s = mod(t, 60);
    m = floor(mod(t-s,3600)/60);
    h = floor((t-s-m*60)/3600);
    outtime = sprintf('%d h %d min %.2f s', h, m, s);
end



