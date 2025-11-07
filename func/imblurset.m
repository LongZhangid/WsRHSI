function u = imblurset(H,x)

[~, ~, num] = size(x);
imset = zeros(size(x));
for jj = 1:num
   
    imout = 0;
    for ii = 1:num
        tmpim = (x(:,:,ii).*H(:,:,ii, jj));
        imout = imout + tmpim;
    end
    imset(:,:,jj) = imout;% + mstd*randn(mm,nn);    
end
u = imset;
end
