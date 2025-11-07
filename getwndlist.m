function wndlist = getwndlist(num,ww,ii)



iia = ii - 1 - floor((ww-1)/2);
iib = ii -1 + ceil((ww-1)/2);

if iia < 0
    iib = iib - iia;
    iia = 0;
end

if iib >= num-1
    iia = iia-(iib-num+1);
    iib = num-1;
end

wndlist = (iia:iib)+1;

