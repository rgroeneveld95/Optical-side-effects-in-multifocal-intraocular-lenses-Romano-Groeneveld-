close all
clc

I1=importdata('intensity05m.mat');
I2=importdata('intensity1m.mat');
I3=importdata('intensity5m.mat');
I4=importdata('intensity10m.mat');
I5=importdata('intensity20m.mat');
rvdraw=importdata('rvdraw.mat');

figure(1)
hold on

semilogy(rvdraw*10^6,I5,'LineWidth',2)
semilogy(rvdraw*10^6,I4,'LineWidth',2)
semilogy(rvdraw*10^6,I3,'LineWidth',2)
semilogy(rvdraw*10^6,I2,'LineWidth',3)
semilogy(rvdraw*10^6,I1,'LineWidth',2)

xlabel('\fontsize{25} r [\mum]');ylabel('\fontsize{25} Intensity \times z_{o}^{2}');
legend('\fontsize{25} z_{o}=20m','\fontsize{25} z_{o}=10m','\fontsize{25} z_{o}=5m','\fontsize{25} z_{o}=1m','\fontsize{25} z_{o}=0.5m','Location','NorthEast')
hold off
set(gca,'FontSize',25)
grid on
ylim([0 4e18]);
xlim([-60 60]);