clc
% rng(2024)
% load simulatedData31-60_std5
% load simulatedData44-50_0
% load simulatedDataTop31-60_0
% load simulatedDataTop44-50_0
% load simulatedData1-31_Pupil_3

% load simulatedData1-31_
% load tmpData
% load simulatedData1-10_
% load simulatedData1-31_Pupil_3_bb64
imsetorg = imset;
yset = imset + 0.001*max(imset,[],'all')*randn(size(imset));  % 0.001

%initialization
sizeF = size(yset(:,:,1));
num = length(f);



eigsDtD = abs(psf2otf([1,-1],sizeF)).^2 + abs(psf2otf([1,-1]',sizeF)).^2;

Hpow = Hset;
Hpow = abs(Hset).^2;


alpha = zeros(num,1);

%加入噪声
% imset = imset + 300*randn(size(imset));   

% flag_alpha = 0;
% 
% if flag_alpha
%     for ii = 1:num
%         tmp = yset(:,:,ii);
%         yset(:,:,ii) = tmp - mean(tmp(:));
%         alpha(ii) = sum(sum(MSsub(:,:,ii)));
%     end
%     alpha = alpha./sum(alpha);
% end


I = yset./num;

psnrI = snrv(I,MSsub)

fftI = I;

% sumHI = zeros(size(imset));
fftim = zeros(size(yset));
for ii = 1:num
    fftI(:,:,ii) = fft2(I(:,:,ii));
    fftim(:,:,ii) = fft2(yset(:,:,ii));     
end





% fftIbuf = fftI;
K = 2;
stdlist=[];
tmp = I(:,:,5) - MSsub(:,:,5);
stdlist = [stdlist std(tmp(:))];

mu1 = 0.05;
mu2 = 0.05;

tic;

for k = 1:K
    k
    
%     fftI = fftIbuf;

    iilist = randperm(num);

%     iilist = 1:num;
    
    for ii = iilist %ii = [1:num]% num-1:1]
%         if k<3
% %             idx = ii;
%             idx = 1:num;
%         else
% %             idx = 1:num;
%             idx = ii;
%         end
        
%         idx = 1:num;
               
        idx = getwndlist(num,1,ii);        
        
        if ii == num
            idneighbors = [ii-2,ii-1];
        elseif ii == 1
            idneighbors = [ii+1,ii+2];
        else
            idneighbors = [ii-1,ii+1];
        end
            
%         end
%         idx = 1:num;
%         idx = ii;
        
        fenmu = 0;
        fenzi = 0;
        
        for ee = idx
            fenmu = fenmu + Hpow(:,:,ii,ee);  %abs(Hset(:,:,ii,ee)).^2;
            tmp = 0;
            
            for jj = 1:num
                if jj ~= ii                    
                    tmp = tmp + Hset(:,:,jj,ee).*fftI(:,:,jj); % - 0.1.*fftI(:,:,ee); 
                end
            end            
            
            tmp = conj(Hset(:,:,ii,ee)).*tmp;

            fenzi = fenzi + conj(Hset(:,:,ii,ee)).*fftim(:,:,ee) - tmp + mu2*(fftI(:,:,idneighbors(1)) + fftI(:,:,idneighbors(2)));             
        end
        
%         fenmu = fenmu + 5*eigsDtD + 0.002;        
        fenmu = fenmu +  0.0005*eigsDtD + mu1+mu2;   %0.0005*eigsDtD     0.005
        
        aa = 0.2;    %0.2   0.15
        fftI(:,:,ii) =(1-aa)*fftI(:,:,ii) +  aa*fenzi./fenmu;
%         fftIbuf(:,:,ii) = fenzi./fenmu;
        



    end
%     fftIbackup = fftI;
    
    tmpI = real(ifft2(fftI(:,:,5)));
    
    tmp = tmpI - (MSsub(:,:,5)-mean(mean(MSsub(:,:,5))));
    stdlist = [stdlist std(tmp(:))];
%     stdlist = [stdlist snr(I,MSsub)];

    

%     fftI = filter(ones(5,1)/5,1,fftI,[],3);
end

% allpower = sum(imsetorg(:))/num;

for ii = 1:num
    I(:,:,ii) = real(ifft2(fftI(:,:,ii)));   
end

time1 = toc;
% for ii = 1:num    
%     if flag_alpha
%         I(:,:,ii) = I(:,:,ii) + allpower*alpha(ii)./length(tmp(:));    
%     end
% end


figure,plot(stdlist)



%%
% 显示rgb彩色图像

[mm,nn,kk] = size(MSsub);

deltL = round((Lb-La+1)/3);

Brance = Flen(La:La+deltL-1,2)*1e-9;
Grance = Flen(La+deltL:La+2*deltL-1,2)*1e-9;
Rrance = Flen(La+2*deltL:Lb,2)*1e-9;

[X,Y,Z] = meshgrid(1:nn,1:mm,lambda);

[Xq,Yq,Zq] = meshgrid(1:nn,1:mm,Rrance);
VR = interp3(X,Y,Z,MSsub,Xq,Yq,Zq);
VRI = interp3(X,Y,Z,I,Xq,Yq,Zq);
VRim = interp3(X,Y,Z,imset,Xq,Yq,Zq);

[Xq,Yq,Zq] = meshgrid(1:nn,1:mm,Grance);
VG = interp3(X,Y,Z,MSsub,Xq,Yq,Zq);
VGI = interp3(X,Y,Z,I,Xq,Yq,Zq);
VGim = interp3(X,Y,Z,imset,Xq,Yq,Zq);

[Xq,Yq,Zq] = meshgrid(1:nn,1:mm,Brance);
VB = interp3(X,Y,Z,MSsub,Xq,Yq,Zq);
VBI = interp3(X,Y,Z,I,Xq,Yq,Zq);
VBim = interp3(X,Y,Z,imset,Xq,Yq,Zq);

Vrgb = zeros(mm,nn,3);
Vrgb(:,:,1) = sum(VR,3);%*1;
Vrgb(:,:,2) = sum(VG,3);%*4.5907;
Vrgb(:,:,3) = sum(VB,3);%*0.0601;
Vrgb = Vrgb/max(Vrgb(:));
figure,imshow(Vrgb,[]);
VrgbI = zeros(mm,nn,3);
VrgbI(:,:,1) = sum(VRI,3);%*1;
VrgbI(:,:,2) = sum(VGI,3);%*4.5907;
VrgbI(:,:,3) = sum(VBI,3);%*0.0601;
VrgbI = VrgbI/max(VrgbI(:));
figure,imshow(VrgbI,[]);

Vrgbim = zeros(mm,nn,3);
Vrgbim(:,:,1) = sum(VRim,3);%*1;
Vrgbim(:,:,2) = sum(VGim,3);%*4.5907;
Vrgbim(:,:,3) = sum(VBim,3);%*0.0601;
Vrgbim = Vrgbim/max(Vrgbim(:));
figure,imshow(Vrgbim,[]);

%%


% clear Brance Rrance Grance VR VRI VRim VG VGI VGim VB VBI VBim


% 
% snr(I,MSsub)
plotSpecLine(MSsub,I,178,306,5);
plotSpecLine(MSsub,I,213,327,5);
plotSpecLine(MSsub,I,101,101,412);

psnro = snrv(I,MSsub)

%% 
% =========================================================================
% Run the deblurring algorithm
% =========================================================================

rng(0)  % random seed, for reproducibility

n_iters = 10;   % 10   200   % number of iterations to solve the denoising subproblem
penalty = @(x) normTVi(x);   % isotropic TV norm as the penalty function
prox_op = @(x,gamma) proxTVi(x,gamma,n_iters);       % proximity operator

x_init = I;  % yset./num; % I; %imgblurT(Hset,fftim); %rand(size(yset));      % random initialization
max_iters = 1500;               % number of iterations
min_iters = 5;
lambda1 = 1e-6; % nature scene: 1e-5  spatial dominate  5e-5    5e-8  TV(xiongan_10) = 6.8532e+07    TV(scene8)=1.3775e+05
y = imgcrop(imset,kernal_size);
y = y + 0.001*max(y,[],'all')*randn(size(y));
% H = zeros(size(fftim));
% xset = zeros(mm, nn, num);
% x = zeros(mm, nn);
% for k = 1:num
%     k
%     for ii = 1:mm
%         for jj = 1:nn
%             
%         end
%     end
%     x = yset(:,:,k);
A  = @(x) imgcrop(real(ifftimg(squeeze(sum(Hset.*fftimg(x), 3)))),kernal_size);    % forward linear operator A
AT = @(x) real(ifftimg(squeeze(sum(conj(permute(Hset, [1, 2, 4, 3])).*fftimg(zeropad(x,kernal_size)), 3))));    % transpose (adjoint operator) of A

[x,n_iters,J_vals,snr_value,~] = TwIST(y,A,lambda1,...   % try ISTA, TwIST, or FISTA
'AT',           AT,...
'initializer',  x_init,...
'prox_op',      prox_op,...
'penalty',      penalty,...
'eta',          2,...
'Lip',          100,... 
'max_iter',     max_iters,...
'min_iter',     min_iters,...
'verbose',      true,...
'objective',    MSsub,...
'psnrori', psnro);
%     I(:,:,k) = x;

% 'objective',    MSsub,...
% 'psnrori', psnro);
% %     I(:,:,k) = x;

% end
[mm,nn,kk] = size(x(:,:,:,1));
I1 = x(:,:,:,1);

% I1 = zeros(mm,nn,kk);
% for ii = 1:num
%     I1(:,:,ii) = real(ifft2(fftI1(:,:,ii)));   
% end


figure,plot(snr_value)



%%
% 显示rgb彩色图像


deltL = round((Lb-La+1)/3);

Brance = Flen(La:La+deltL-1,2)*1e-9;
Grance = Flen(La+deltL:La+2*deltL-1,2)*1e-9;
Rrance = Flen(La+2*deltL:Lb,2)*1e-9;

[X,Y,Z] = meshgrid(1:nn,1:mm,lambda);

[Xq,Yq,Zq] = meshgrid(1:nn,1:mm,Rrance);
% VR = interp3(X,Y,Z,MSsub,Xq,Yq,Zq);I1
VRI = interp3(X,Y,Z,I1,Xq,Yq,Zq);
% VRim = interp3(X,Y,Z,imset,Xq,Yq,Zq);

[Xq,Yq,Zq] = meshgrid(1:nn,1:mm,Grance);
% VG = interp3(X,Y,Z,MSsub,Xq,Yq,Zq);
VGI = interp3(X,Y,Z,I1,Xq,Yq,Zq);
% VGim = interp3(X,Y,Z,imset,Xq,Yq,Zq);

[Xq,Yq,Zq] = meshgrid(1:nn,1:mm,Brance);
% VB = interp3(X,Y,Z,MSsub,Xq,Yq,Zq);
VBI = interp3(X,Y,Z,I1,Xq,Yq,Zq);
% VBim = interp3(X,Y,Z,imset,Xq,Yq,Zq);

% Vrgb = zeros(mm,nn,3);
% Vrgb(:,:,1) = sum(VR,3);%*1;
% Vrgb(:,:,2) = sum(VG,3);%*4.5907;
% Vrgb(:,:,3) = sum(VB,3);%*0.0601;
% Vrgb = Vrgb/max(Vrgb(:));
% figure,imshow(Vrgb,[]);
VrgbI = zeros(mm,nn,3);
VrgbI(:,:,1) = sum(VRI,3);%*1;
VrgbI(:,:,2) = sum(VGI,3);%*4.5907;
VrgbI(:,:,3) = sum(VBI,3);%*0.0601;
VrgbI = VrgbI/max(VrgbI(:));
figure,imshow(VrgbI,[]);

% Vrgbim = zeros(mm,nn,3);
% Vrgbim(:,:,1) = sum(VRim,3);%*1;
% Vrgbim(:,:,2) = sum(VGim,3);%*4.5907;
% Vrgbim(:,:,3) = sum(VBim,3);%*0.0601;
% Vrgbim = Vrgbim/max(Vrgbim(:));
% figure,imshow(Vrgbim,[]);

%%





% 
% snr(I,MSsub)
plotSpecLine(MSsub,I1,178,306,5);
plotSpecLine(MSsub,I1,213,327,5);
plotSpecLine(MSsub,I1,101,101,412);
%% 
% =========================================================================
% Auxiliary functions
% =========================================================================

function u = imgcrop(x,cropsize)
% =========================================================================
% Crop the central part of the image.
% -------------------------------------------------------------------------
% Input:    - x        : Original image.
%           - cropsize : Cropping pixel number along each dimension.
% Output:   - u        : Cropped image.
% =========================================================================
u = x(cropsize+1:end-cropsize,cropsize+1:end-cropsize,:);
end

function u = zeropad(x,padsize)
% =========================================================================
% Zero-pad the image.
% -------------------------------------------------------------------------
% Input:    - x        : Original image.
%           - padsize  : Padding pixel number along each dimension.
% Output:   - u        : Zero-padded image.
% =========================================================================
u = padarray(x,[padsize,padsize],0);
end

function u = fftimg(x)
% =========================================================================
% Zero-pad the image.
% -------------------------------------------------------------------------
% Input:    - x        : Original image.
%           - padsize  : Padding pixel number along each dimension.
% Output:   - u        : Zero-padded image.
% =========================================================================
[m, n, num] = size(x);
u = zeros(m, n, num);
for j = 1:num
    u(:,:,j) = fft2(x(:,:,j));
end
end

function u = ifftimg(x)
% =========================================================================
% Zero-pad the image.
% -------------------------------------------------------------------------
% Input:    - x        : Original image.
%           - padsize  : Padding pixel number along each dimension.
% Output:   - u        : Zero-padded image.
% =========================================================================
[m, n, num] = size(x);
u = zeros(m, n, num);
for j = 1:num
    u(:,:,j) = ifft2(x(:,:,j));
end
end