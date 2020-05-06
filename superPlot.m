%% plot AMS Agglomeration - Fivespot Structured
% pick a k 1, 1.5 2 3 4 5
k = 3;
auxmat = load('ams-agglomaration-fivespot-structured-cf.txt');

type1 = auxmat(:,2) == 1;
type2 = ~type1;
kref = auxmat(:,3) == k;
figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 500 800])
set(figHandle, 'name', ['Production Curves'],'NumberTitle','off');
whitebg('white')

estrutura = 'Agglomeration - Structured Grid ';

subplot(2,1,1)
% normal L2
minCV = min(auxmat(:,1));
maxCV = max(auxmat(:,1));

plot(auxmat(type1 & kref,1), auxmat(type1 & kref,4), 'DisplayName','Old Dual', 'color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2 & kref,1), auxmat(type2 & kref,4), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show','Location','southeast')
title(['Five Spot - ', estrutura,  '- Kyy/Kxx = ', num2str(k)],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('L2 Norm (%)','FontSize',10,'FontWeight','normal','Color',[0 0 0])
% normal Linf
subplot(2,1,2)

plot(auxmat(type1 & kref,1), auxmat(type1 & kref,5), 'DisplayName','Old Dual','color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2 & kref,1), auxmat(type2 & kref,5), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show')
title(['Five Spot - ', estrutura,  '- Kyy/Kxx = ', num2str(k)],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('LInfinity Norm (%)','FontSize', 10,'FontWeight','normal','Color',[0 0 0])

saveas(figHandle,'test.png')
%% plot AMS Agglomeration - Fivespot Unstructured
% pick a k 1, 1.5 2 3 4 5
for k = [1 1.5 2 3 4 5]
%k = 1;

auxmat = load('ams-agglomaration-fivespot-unstructured-cf.txt');

type1 = auxmat(:,2) == 1;
type2 = ~type1;
kref = auxmat(:,3) == k;
figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 500 800])
set(figHandle, 'name', ['Production Curves'],'NumberTitle','off');
whitebg('white')
estrutura = 'Agglomeration - Unstructured Grid ';
subplot(2,1,1)
% normal L2
minCV = min(auxmat(:,1));
maxCV = max(auxmat(:,1));
plot(auxmat(type1 & kref,1), auxmat(type1 & kref,4), 'DisplayName','Old Dual', 'color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2 & kref,1), auxmat(type2 & kref,4), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show','Location','southeast')
title(['Five Spot - ', estrutura,  '- Kyy/Kxx = ', num2str(k)],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('L2 Norm (%)','FontSize',10,'FontWeight','normal','Color',[0 0 0])
% normal Linf
subplot(2,1,2)
plot(auxmat(type1 & kref,1), auxmat(type1 & kref,5), 'DisplayName','Old Dual','color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2 & kref,1), auxmat(type2 & kref,5), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show')
title(['Five Spot - ', estrutura,  '- Kyy/Kxx = ', num2str(k)],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('LInfinity Norm (%)','FontSize', 10,'FontWeight','normal','Color',[0 0 0])

gname = ['uns-agg-' num2str(k) '.png']
saveas(gcf,gname )
end

%% plot AMS Hex - Linear Structured
% pick a k 1, 1.5 2 3 4 5
k = 1;
for k = [1 1.5 2 3 4 5]
auxmat = load('amshexlinear.txt');
auxvec = load('dupar.txt');
%auxmat(:,1) = auxvec;
type1 = auxmat(:,2) == 1;
type2 = ~type1;
kref = auxmat(:,3) == k;
figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 500 800])
set(figHandle, 'name', ['Production Curves'],'NumberTitle','off');
whitebg('white')
estrutura = 'Honeycomb - Unstructured Grid ';
subplot(2,1,1)
% normal L2
minCV = min(auxmat(:,1));
maxCV = max(auxmat(:,1));
minCV = min(auxvec);
maxCV = max(auxvec);
plot(auxvec, auxmat(type1 & kref,4), 'DisplayName','Old Dual', 'color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxvec, auxmat(type2 & kref,4), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show','Location','southeast')
title(['Flow Channel - ', estrutura,  '- Kyy/Kxx = ', num2str(k)],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('L2 Norm (%)','FontSize',10,'FontWeight','normal','Color',[0 0 0])
% normal Linf
subplot(2,1,2)
plot(auxvec, auxmat(type1 & kref,5), 'DisplayName','Old Dual','color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxvec, auxmat(type2 & kref,5), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show')
title(['Flow Channel - ', estrutura,  '- Kyy/Kxx = ', num2str(k)],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])  
ylabel('LInfinity Norm (%)','FontSize', 10,'FontWeight','normal','Color',[0 0 0])

gname = ['uns-hex-' num2str(k) '.png']
saveas(gcf,gname )
end



%% plot AMS Hex - Five Spot - 01
% pick a k 1, 1.5 2 3 4 5
k = 1;
for k = 1
auxmat = load('ams-hex-fivespot-unstructured-hh01.txt');
%auxmat(:,1) = auxvec;
type1 = auxmat(:,2) == 1;
type2 = ~type1;
kref = type1& type2;%auxmat(:,3) == k;
figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 500 800])
set(figHandle, 'name', ['Production Curves'],'NumberTitle','off');
whitebg('white')
estrutura = 'Honeycomb - Unstructured Grid ';
subplot(2,1,1)
% normal L2
minCV = min(auxmat(:,1));
maxCV = max(auxmat(:,1));

plot(auxmat(type1,1), auxmat(type1,3), 'DisplayName','Old Dual', 'color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,3), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show','Location','southeast')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('L2 Norm (%)','FontSize',10,'FontWeight','normal','Color',[0 0 0])
% normal Linf
subplot(2,1,2)
plot(auxmat(type1,1), auxmat(type1,4), 'DisplayName','Old Dual','color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,4), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])  
ylabel('LInfinity Norm (%)','FontSize', 10,'FontWeight','normal','Color',[0 0 0])

gname = ['uns-hex-symmetric-heterogeneous-reservoir.png'];
saveas(gcf,gname )
end


%% plot AMS Hex - Five Spot - 01
% pick a k 1, 1.5 2 3 4 5
k = 1;
for k = 1
auxmat = load('ams-agglomeration-fivespot-structured-hh01.txt');
%auxmat(:,1) = auxvec;
type1 = auxmat(:,2) == 1;
type2 = ~type1;
kref = type1& type2;%auxmat(:,3) == k;
figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 500 800])
set(figHandle, 'name', ['Production Curves'],'NumberTitle','off');
whitebg('white')
estrutura = 'Agglomeration- Structured Grid ';
subplot(2,1,1)
% normal L2
minCV = min(auxmat(:,1));
maxCV = max(auxmat(:,1));

plot(auxmat(type1,1), auxmat(type1,3), 'DisplayName','Old Dual', 'color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,3), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show','Location','southeast')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('L2 Norm (%)','FontSize',10,'FontWeight','normal','Color',[0 0 0])
% normal Linf
subplot(2,1,2)
plot(auxmat(type1,1), auxmat(type1,4), 'DisplayName','Old Dual','color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,4), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])  
ylabel('LInfinity Norm (%)','FontSize', 10,'FontWeight','normal','Color',[0 0 0])

gname = ['uns-agg-symmetric-heterogeneous-reservoir.png'];
saveas(gcf,gname )
end



%% plot AMS Hex - Five Spot - 02
% pick a k 1, 1.5 2 3 4 5
k = 1;
for k = 1
auxmat = load('ams-agglomeration-fivespot-structured-hh02.txt');
%auxmat(:,1) = auxvec;
type1 = auxmat(:,2) == 1;
type2 = ~type1;
kref = type1& type2;%auxmat(:,3) == k;
figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 500 800])
set(figHandle, 'name', ['Production Curves'],'NumberTitle','off');
whitebg('white')
estrutura = 'Agglomeration- Structured Grid ';
subplot(2,1,1)
% normal L2
minCV = min(auxmat(:,1));
maxCV = max(auxmat(:,1));

plot(auxmat(type1,1), auxmat(type1,3), 'DisplayName','Old Dual', 'color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,3), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show','Location','southeast')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('L2 Norm (%)','FontSize',10,'FontWeight','normal','Color',[0 0 0])
% normal Linf
subplot(2,1,2)
plot(auxmat(type1,1), auxmat(type1,4), 'DisplayName','Old Dual','color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,4), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])  
ylabel('LInfinity Norm (%)','FontSize', 10,'FontWeight','normal','Color',[0 0 0])

gname = ['struc-agg-aleatory-heterogeneous-reservoir.png'];
saveas(gcf,gname )
end


%% plot AMS Hex - Five Spot - 02
% pick a k 1, 1.5 2 3 4 5
k = 1;
for k = 1
auxmat = load('ams-agglomeration-fivespot-unstructured-hh02.txt');
%auxmat(:,1) = auxvec;
type1 = auxmat(:,2) == 1;
type2 = ~type1;
kref = type1& type2;%auxmat(:,3) == k;
figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 500 800])
set(figHandle, 'name', ['Production Curves'],'NumberTitle','off');
whitebg('white')
estrutura = 'Agglomeration- Structured Grid ';
subplot(2,1,1)
% normal L2
minCV = min(auxmat(:,1));
maxCV = max(auxmat(:,1));

plot(auxmat(type1,1), auxmat(type1,3), 'DisplayName','Old Dual', 'color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,3), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show','Location','southeast')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('L2 Norm (%)','FontSize',10,'FontWeight','normal','Color',[0 0 0])
% normal Linf
subplot(2,1,2)
plot(auxmat(type1,1), auxmat(type1,4), 'DisplayName','Old Dual','color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,4), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])  
ylabel('LInfinity Norm (%)','FontSize', 10,'FontWeight','normal','Color',[0 0 0])

gname = ['unstruc-agg-aleatory-heterogeneous-reservoir.png'];
saveas(gcf,gname )
end



%% plot AMS Hex - Five Spot - 02
% pick a k 1, 1.5 2 3 4 5
k = 1;
for k = 1
auxmat = load('ams-hex-fivespot-unstructured-hh02.txt');
%auxmat(:,1) = auxvec;
type1 = auxmat(:,2) == 1;
type2 = ~type1;
kref = type1& type2;%auxmat(:,3) == k;
figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 500 800])
set(figHandle, 'name', ['Production Curves'],'NumberTitle','off');
whitebg('white')
estrutura = 'Honeycomb - Unstructured Grid ';
subplot(2,1,1)
% normal L2
minCV = min(auxmat(:,1));
maxCV = max(auxmat(:,1));

plot(auxmat(type1,1), auxmat(type1,3), 'DisplayName','Old Dual', 'color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,3), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show','Location','southeast')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])
ylabel('L2 Norm (%)','FontSize',10,'FontWeight','normal','Color',[0 0 0])
% normal Linf
subplot(2,1,2)
plot(auxmat(type1,1), auxmat(type1,4), 'DisplayName','Old Dual','color',[0 0.2969 0.5977],'LineWidth',2.5);
hold on
plot(auxmat(type2,1), auxmat(type2,4), 'DisplayName','New Dual', 'color',[0.782 0 0.3],'LineWidth',2.5);
xlim([minCV maxCV ])
legend('show')
title(['1/4 of Five Spot - ', estrutura],'FontSize', 10)
xlabel('Number of Coarse Volumes','FontSize',10,'FontWeight','normal','Color',[0 0 0])  
ylabel('LInfinity Norm (%)','FontSize', 10,'FontWeight','normal','Color',[0 0 0])

gname = ['unstruc-hex-aleatory-heterogeneous-reservoir.png'];
saveas(gcf,gname )
end