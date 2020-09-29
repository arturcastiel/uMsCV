figHandle = figure;
figure(figHandle);
% 
%en = 1000;
set(figHandle, 'Position', [0 0 800 800])
set(figHandle,'color','w')
%axis off
% hold on
control = importdata(strcat(superFolder,'\control.dat'));
residos = importdata(strcat(superFolder,'\res.dat'));

graficos = cell(size(control,1),2);

st = 1;
for ii = 1:size(control,1)
   it = control(ii,1);
   num = control(ii,2);
   yy = residos(st: st+num-1);
   xx = linspace(0,it, num);
   ref = ismember(xx,0:num);
   graficos{ii,1} = xx(ref);
   graficos{ii,2} = yy(ref);    
   st = st + num;
   semilogy(xx(ref),yy(ref), 'LineWidth', 3);
   hold on
   grid on

    %%semilogy(xx,yy,'LineWidth', 3);
    xlim([0 max(xx(ref))])
    legend('BiCGSTAB','BiCGSTAB-with-osc')
end