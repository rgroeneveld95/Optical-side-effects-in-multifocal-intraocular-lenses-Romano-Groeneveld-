close all
clc

I1=importdata('axialintensity05v2.mat');
I2=importdata('axialintensity10v2.mat');
I3=importdata('axialintensity15v2.mat');
I4=importdata('axialintensity20v2.mat');
I5=importdata('axialintensity25v2.mat');
I6=importdata('axialintensity30v2.mat');
zcor=importdata('axialcoordinatev2.mat');

figure(1)
hold on
plot(zcor*1000,I6,'LineWidth',1.5);
plot(zcor*1000,I5,'LineWidth',1.5);
plot(zcor*1000,I4,'LineWidth',1.5);
plot(zcor*1000,I3,'LineWidth',1.5);
plot(zcor*1000,I2,'LineWidth',1.5);
plot(zcor*1000,I1,'LineWidth',1.5);
xlabel('\fontsize{25} z_i [mm]'); 
ylabel('\fontsize{25} Axial intensity / a_{e}^{2}');
leg = legend('show');
title(leg,'\fontsize{25} pupil diameter')
legend('\fontsize{25} 6 mm','\fontsize{25} 5 mm','\fontsize{25} 4mm','\fontsize{25} 3mm','\fontsize{25} 2mm','\fontsize{25} 1mm')
hold off
% xlim([18 22]);
grid on
set(gca,'FontSize',25)