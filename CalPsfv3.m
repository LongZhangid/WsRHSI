function psf = CalPsfv3(D,f,M,N,lambda,z,pixs,Pimg)
% Calculate the point spread function. 
% D is the Caliber of optical system.
% f is the Focal length of the optical system.It can be a array.
% M is the height of the psf.
% N is the width of the psf.
% lambda is the wave length. It can be a array.
% z is the distance from the exit pupil to the ccd. It can be a array.
% pixs is the pixel size

NN = max([512, M, N]);
num = length(lambda);

% make the pupil matrix.
w = NN;
W = zeros(w,w,num);
psf= zeros(M,N,num);
P = zeros(w,w);

R = f/D.*lambda./pixs;

dxx = 2.0./(w-1);
	xx = -1;
	yy = 1;
for (c0=0:1:w-1)
    for (c1=0:1:w-1)
        ro = sqrt(xx.*xx + yy.*yy);	
        W(c0+1,c1+1,:) = pi.*(z-f).*((ro*D./f).^2)./lambda./4;
        xx=xx+dxx;     
%         P(c0+1,c1+1,ro>1)=0;
%         W(c0+1,c1+1,ro>1)=0;
%         P(c0+1,c1+1,ro<=1)=1;  
        if ro>1
%             P(c0+1,c1+1)=0;
            W(c0+1,c1+1,:)=0;
        else
%             P(c0+1,c1+1)=1;
        end  
    end
    xx = -1;
    yy=yy-dxx;
end
P = imresize(Pimg,[NN,NN]);

for ii = 1:num
   
    Gpupil=P.*exp(i*(W(:,:,ii)));
    Gpupil=fftshift(Gpupil);
    hh = ifft2(Gpupil);
    tmppsf=fftshift(hh.*conj(hh));
%     tmppsf = tmppsf./sum(sum(tmppsf));
    tmppsf = imresize(tmppsf,R(ii)/1);
    [w1,w2] = size(tmppsf);
    if w>w1
        tmppsf = padarray(tmppsf,ceil([(w-w1)/2,(w-w2)/2]));
        ww = w;
    else
        ww = w1;
    end
    n1 = floor((ww-N)/2) + 1;
    n2 = n1 + N - 1;
    m1 = floor((ww-M)/2) + 1;
    m2 = m1 + M - 1;
    tt = tmppsf(m1:m2,n1:n2);    
    psf(:,:,ii) = tt./sum(tt(:));
end



