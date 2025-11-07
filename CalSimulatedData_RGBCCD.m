clc
clear
close all
%Psfset(:,:,ii,jj)
%获取的第jj个图像，的第ii个波段的点扩散函擿

load scene10  % ref_mosteiro4bb_reg1_lax  %galaxy   % ref_cbrufefields1bb_reg1_lax  % ref_cyflower1bb_reg1     % ref_ruivaes1bb_reg1_lax    % ref_crown3bb_reg1_lax    % Houston  % xiongan    % scene5  % 0069  % 0020   scene1   scene3  scene8

pixel_size = 12*1e-6;

MS = double(reflectances(211:600,211:600,:))/10000;  % reflectances(211:600,211:600,:)  %reflectances(161:672,161:672,:);   XiongAn(300:600,300:600,1:144)    Houston(:,200:548,1:144)  hypercube(149:620,277:748,1:31)
% MS = MS ./ max(MS,[], 'all');
% figure,imshow(MS)

kernal_size = 100;
tmp_size = size(MS,1,2) + 2*kernal_size;
tmp = zeros([tmp_size, size(MS,3)]);
% tmp(498:574,484:563,:) = MS(498:574,484:563,:);
% tmp(298:574,284:563,:) = MS(298:574,284:563,:);
tmp(kernal_size+1:end-kernal_size,kernal_size+1:end-kernal_size,:) = MS;  %(kernal_size+1:end-kernal_size,kernal_size+1:end-kernal_size,:);

% MS = tmp;

lmd = 400:10:700;   %400:10:700;   500:1:530;
ff = 60000./lmd;
Flen = [(1:31)',lmd',ff'];

% ii1 = 800;
% ii2 = ii1+256-1;

% ii1 = 1;
% ii2 = 512;

% MS = MS(ii1:ii2,1:512,:);

[mm,nn,~] = size(tmp);

Pimg = double(imread('Pupil.png'));
Pimg = Pimg(:,:,1)/255;
tmph = fspecial('gaussian',200,60);
Pimg = imfilter(Pimg,tmph,'conv');


%Flen 第一列为序号；第二列为波长，单位nm；第三列为焦距，单位mm

% La = 1;
La = 11;
Lb = 20;
% La = 6;
% Lb = 25;
mstd = 0;
MSsub = double(tmp(:,:,La:Lb));
lambda = Flen(La:Lb,2)*1e-9;
f = Flen(La:Lb,3)*1e-3;
num = length(lambda);

pw = min(mm,nn);
Psfset = zeros(pw,pw,num,num);
imset = zeros(mm,nn,num);

% imset15 = zeros(size(MSsub));

for jj = 1:num
    jj
    z = f(jj)*ones(num,1);
%     psf = CalPsf(0.01,f,256,256,lambda,z);
%     psf = CalPsfv2(0.01,f,pw,pw,lambda,z,pixel_size);
    psf = CalPsfv3(0.005/0.7,f,pw,pw,lambda,z,pixel_size,Pimg);  % D: 0.01/0.7
    Psfset(:,:,:,jj) = psf;
    imout = 0;
    for ii = 1:length(lambda)
%         imout = imout + imfilter(MSsub(:,:,ii),psf(:,:,ii),'circular','conv');        
%         imout = imout + ifft2(fft2(MSsub(:,:,ii)).*fft2(fftshift(psf(:,:,ii))));
        tmpH = psf2otf(psf(:,:,ii),[mm,nn]);
        tmpim = ifft2(fft2(MSsub(:,:,ii)).*tmpH);
%         if jj == 15
%             imset15(:,:,ii) = tmpim;
%         end        
        imout = imout + tmpim;
    end
    imset(:,:,jj) = imout;% + mstd*randn(mm,nn);    
end


Hset = zeros([mm,nn,num,num]);
for ii = 1:num
    for jj = 1:num
        Hset(:,:,ii,jj) = psf2otf(Psfset(:,:,ii,jj),[mm,nn]);
    end
end

% save(['simulatedData',num2str(La),'-',num2str(Lb),'_Pupil_3','.mat'],'imset','Hset','MSsub','f','lambda','-v7.3')


clear MS imout tmpim tmpH tmp reflectances lmd Pimg psf tmph z XiongAn Houston hypercube 0020 0031 0028 0069 0039

% save('tmpData.mat','imset','Psfset','MSsub','f','lambda')




