function plotSpecLine1(A,B,x,y,w)
% 显示光谱立方体 A 和 B 在x和y位置的窗口宽为w的光谱曲线
[n,m,k] = size(A);
a = A(y:y+w-1,x:x+w-1,:);
a = reshape(a,[w*w,k]);
b = B(y:y+w-1,x:x+w-1,:);
b = reshape(b,[w*w,k]);
if w>1
a = mean(a);
b = mean(b);
end
band = linspace(621, 630, k);
figure,plot(band,a);
hold on, plot(band,b);
% 设置x轴和y轴
xlim([621 630]); % x轴从500到700
ylim([0 1.5]); % y轴从0到1
% 设置x轴和y轴刻度间隔为
xticks(621:1:630);
yticks(0:0.2:1.5);
legend('Original data','Reconstructed data', 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Wavelength(nm)', 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Intensity', 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
std(a-b)
% b = filter(ones(5,1)/5,1,b);
% figure,plot(a);
% hold on, plot(b);