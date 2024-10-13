%%%%the original scale picture
clear all;
close all;
FS=18;
%dat=load('Wave_Group6.dat');
dat=load('DNA.dat');
figure(1);
set(gcf, 'position', [100 100 600 500],'color','w');

loglog(dat(2:end,2)/60,dat(2:end,3),'*','markersize',8,'color','b')
hold on
%errorbar(dat(2:end,2),dat(2:end,3),dat(2:end,4))
%loglog(dat(2:end,2),dat(2:end,3),'o','markersize',11)
%errorbar(dat(2:end,2),dat(2:end,3),dat(2:end,4))
errorbar(dat(2:end,2)/60,dat(2:end,3),dat(2:end,4),'color','g')
hold on

x=linspace(5,35,100);
%y=13*x.^0.35;
y=2.4*x.^0.28;
loglog(x,y,'r--','linewidth',2)
xlabel('time (hour)')
ylabel('spatial scale (mm)')

text(20,4,'$l\sim t^{0.28}$','Interpreter','latex','fontsize',16)

xlabel('Time, $t$ [min]','Interpreter','latex');
ylabel('Spatial scale, $\ell$ [mm]','Interpreter','latex');
%ylim([10 35]);
 ylim([3 8]);
%xlim([0.9 400]);
xlim([5 150]);
xticks([5 15 30  60  150]);
%yticks([ 15  20 25 30 40]);
yticks([3 4 5 6 8]);


set(gca,'XScale','log','YScale','log');
set(gca,'fontsize',FS,'linewidth',2,'xminortick','on','yminortick','on',...
    'ticklength',[0.025 0.01]);
set(gca,'FontName','Times'); set(gcf,'Color',[1,1,1]);

% filename = 'DNA_local_wave_realdata_1118.pdf';
% 
% saveas(gcf, filename);
% filename2 = 'DNA_local_wave_realdata_1118.jpg';
% saveas(gcf, filename2);