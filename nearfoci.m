%finding the near focus

%Run master script
Optical_side_effects_in_multifocal_intraocular_lenses_mep; 
close all
clc

%=============================================================================================
%CALCULATION
%=============================================================================================

k1=(n_a-n_air)/Rc;
k2=(n_iol-n_a)/Ra;
k3=(n_v-n_iol)/Rp;

Mp=[1 -k2; 0 1];
Md=[1 0;d_iol/n_iol 1];
Ma=[1 -k3; 0 1];

%Transfer matrices of the corneal lens and IOL lens
MV1=[1 -k1;0 1];
MV3V2=Mp*Md*Ma;

%Optical powers of the corneal lens and IOL lens
Pc=-MV1(1,2);
P_iol=-MV3V2(1,2);

%T1 & T2
T1=n_a/MV3V2(1,2)*(MV3V2(1,1)-1);
T2=n_v*(MV3V2(2,2)/MV3V2(1,2)*(MV3V2(1,1)-1)-MV3V2(2,1));

D=V2-V1+T1;


s01=linspace(-50,50,1000000);

si1=n_a*s01./(Pc*s01-n_air);
s02=D-si1;
si2=n_v*s02./(P_iol*s02-n_a);
di2=si2+H2-z_e;

s011=geo_image_inverse(si2,n_air,n_a,n_iol,n_v,Rc,Ra,Rp,d_iol,V1,V2);




figure(1)

hold on
plot(s01*100,di2*1000,'.b','LineWidth',2)
plot(s01*100,19.33e-3*ones(1,length(di2))*1000,'r','LineWidth',1)
plot(s01*100,19.6e-3*ones(1,length(di2))*1000,'color','#fc6500','LineWidth',1)
plot(-0.33505*ones(1,length(di2))*100,di2*1000,'r','LineWidth',1)
plot(-0.43375*ones(1,length(di2))*100,di2*1000,'color','#fc6500','LineWidth',1)
text(-0.33505*100,20.5,{['\fontsize{25} z_{o}=',num2str(-33.51),'cm']},'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','color','red');
text(-0.5*100,20.5,{['\fontsize{25} z_{o}=',num2str(-43.38),'cm']},'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','color','#fc6500');
hold off
legend('\fontsize{25} z_{i} vs. z_{o}','\fontsize{25} SN6AD3','\fontsize{25} SN6AD1')
xlabel('\fontsize{25} z_{o} [cm]');ylabel('\fontsize{25} z_{i} [mm]')
ylim([17 22]);
xlim([-2.5*100 2.5*100]);
set(gca,'FontSize',25)
grid on