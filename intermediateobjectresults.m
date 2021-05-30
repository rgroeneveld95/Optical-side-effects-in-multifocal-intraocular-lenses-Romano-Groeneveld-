%Run master script
Optical_side_effects_in_multifocal_intraocular_lenses_mep; 
close all
clc
close all
clc

%=============================================================================================
%GEOMETRICAL OPTICS: DETERMINATION OF IMAGE DISTANCE CALCULATION
%=============================================================================================

%SYMBOL EXPLANATION

%Vertices

%V1 = the vertex of the cornea
%V2 = the vertex of the anterior surface of the smooth IOL lens
%V3 = the vertex of the posterior surface of the smooth IOL lens

%Distances 

%So1 = object distance w.r.t. primary principal plane H1
%T1 = distance of the primary principal plane H1 to the primary vertex V1
%d = anterior chamber width measured from vertex V1 to vertex V2
%d_iol = width of the IOL measured from vertex V2 to vertex V3
%T2 = distance of the secondary principal plane H2 to the last vertex V3
%Si2 = image distance w.r.t. secondary principal plane H2

%Radii of curvature

%Rc = corneal curvature
%Ra = anterior surface curvature of the smooth IOL lens
%Rp = posterior surface curvature of the smooth IOL lens

%Refractive indices

%n_air = refractive index of air
%n_a = refractive index aqueous humor
%n_iol = refractive index of the IOL
%n_v = refractive index of the vitreous humor

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


%Plotting graphs

figure(1)
hold on
plot(s01,si1*1000,'.b','LineWidth',2)
hold off
legend('\fontsize{25} S_{1i} vs. S_{1o}')
xlabel('S_{1o} [m]');ylabel('S_{1i} [mm]')
xlim([-50 50]);
ylim([-120 180]);
set(gca,'FontSize',25)
grid on


figure(2)
hold on
plot(si1*1000,s02*1000,'.b','LineWidth',2)
hold off
legend('\fontsize{25} S_{2o} vs. S_{1i}')
xlabel('\fontsize{25} S_{1i} [mm]');ylabel('\fontsize{25} S_{2o} [mm]')
xlim([-120 180]);
ylim([-200 150]);
set(gca,'FontSize',25)
grid on


figure(3)
hold on
plot(s02*1000,si2,'.b','LineWidth',2)
hold off
legend('\fontsize{25} S_{2i} vs. S_{2o}')
xlabel('\fontsize{25} S_{2o} [mm]');ylabel('\fontsize{25} S_{2i} [m]')
set(gca,'FontSize',25)
grid on
xlim([-50 180]);
ylim([-6 6]);

figure(4)
TF1x = islocalmin(di2);
TF2x = islocalmax(di2);
point1=(s01(TF1x)+s01(TF2x))/2;
point2=(di2(1)+di2(end))/2;

hold on
plot(s01,di2*1000,'.b','LineWidth',2)
text(point1,point2,{['\fontsize{25} f1=',num2str(point1*1000,5),'mm']},'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','color','blue');
text(15,point2*1000,{['\fontsize{25} f2=',num2str(point2*1000,5),'mm']},'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','color','blue');
plot(s01,d_r*ones(1,length(di2))*1000,'--r','LineWidth',2)
hold off
legend('\fontsize{25} z_{i} vs. z_{o}','\fontsize{25} z_{r}')
xlabel('\fontsize{25} z_{o} [m]');ylabel('\fontsize{25} z_{i} [mm]')
ylim([0 40]);
xlim([-30 30]);
set(gca,'FontSize',25)
grid on





