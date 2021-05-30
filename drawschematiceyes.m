%Run master script
Optical_side_effects_in_multifocal_intraocular_lenses_mep; 
close all
clc


%Samples
Ndraw=200;                                  %No. of samples
d=15e-3;                                    %Calculation domain in radial direction
rsample=linspace(-d,d,Ndraw);               %Sample space in radial direction
asample=linspace(0,2*pi,Ndraw);             %Sample space in azimuthal direction
r_c=3e-3;
scale=1000;

%Schematic eyes components

%Surface eye
z_eye=d_eye/2+d_eye/2*cos(asample);
rho_eye=d_eye/2*sin(asample);

%Surface cornea
rho_cornea=linspace(-r_c,r_c,Ndraw);
z_cornea=V1+r_c-sqrt(r_c^2-rho_cornea.^2);

%Smooth anterior and posterior surface IOL
rho_IOL=linspace(-r_iol,r_iol,Ndraw);
z_aIOL=V2+Ra-sqrt(Ra^2-rho_IOL.^2);
z_pIOL=V3+Rp+sqrt(Rp^2-rho_IOL.^2);

%Aperture and exit pupil
rho_apup=linspace(r_pupil,8.6e-3,Ndraw);
rho_apdown=linspace(-8.6e-3,-r_pupil,Ndraw);
rho_exup=rho_apup*Mt;
rho_exdown=rho_apdown*Mt;

%Rays
z_rays=linspace(z_e,d_eye,Ndraw);
r_rays=(z_rays-d_eye)*-(r_e/(d_eye-z_e));
z_imaging=linspace(z_e,20e-3,Ndraw);
r_imaging=(z_imaging-20e-3)*-(r_e/(20e-3-z_e));

%Plotting

%Close up
figure(1)
hold on

%Rays
p1=plot(z_rays*scale,r_rays*scale,'-.k','LineWidth',1.5);
p2=plot(z_rays*scale,-r_rays*scale,'-.k','LineWidth',1.5);

%Vertices
p3=plot(V1*ones(1,Ndraw)*scale,rsample*scale,'--r','LineWidth',1.5);
p4=plot(V2*ones(1,Ndraw)*scale,rsample*scale,'--','color','#b22929','LineWidth',1.5);
p5=plot(V3*ones(1,Ndraw)*scale,rsample*scale,'--','color','#0095c9','LineWidth',1.5);
p6=plot(d_eye*ones(1,Ndraw)*scale,rsample*scale,'--k','LineWidth',1.5);
p17=plot(z_e*ones(1,Ndraw)*scale,rsample*scale,'--b','LineWidth',1.5);


%Principal planes
p7=plot(H1*ones(1,Ndraw)*scale,rsample*scale,'m','LineWidth',1.5);
p8=plot(H2*ones(1,Ndraw)*scale,rsample*scale,'color','#a0dfd4','LineWidth',1.5);

%Aperture & Exit pupil
p9=plot(V2*ones(1,Ndraw)*scale,rho_apup*scale,'color','#FFA500','LineWidth',2);
p10=plot(V2*ones(1,Ndraw)*scale,rho_apdown*scale,'color','#FFA500','LineWidth',2);
p11=plot(z_e*ones(1,Ndraw)*scale,rho_exup*scale,'b','LineWidth',2);
p12=plot(z_e*ones(1,Ndraw)*scale,rho_exdown*scale,'b','LineWidth',2);

%Contours
p13=plot(z_eye*scale,rho_eye*scale,'k','LineWidth',2);
p14=plot(z_cornea*scale,rho_cornea*scale,'r','LineWidth',2);
p15=plot(z_aIOL*scale,rho_IOL*scale,'color','#b22929','LineWidth',2);
p16=plot(z_pIOL*scale,rho_IOL*scale,'color','#0095c9','LineWidth',2);

hold off
xlabel('\fontsize{25} z [mm]');ylabel('\fontsize{25} r [mm]');
legend([p1 p3 p4 p5 p17 p7 p8 p9 p11 p14 p15 p16],{'Imaging cone','V_{c}','V_{a}','V_{p}','z_{e}','H_{1}','H_{2}','AS','EP','LC','AI','PI'},'location','best');
set(gca,'FontSize',25)
grid on
xlim([-0.25 6]);
ylim([-3.5 3.5]);


%Whole eye
figure(2)
hold on

%Rays
p1=plot(z_rays*scale,r_rays*scale,'-.k','LineWidth',1.5);
plot(z_rays*scale,-r_rays*scale,'-.k','LineWidth',1.5);

%Rays arbitrary image
p2=plot(z_imaging*scale,r_imaging*scale,'-.g','LineWidth',1.5);
plot(z_imaging*scale,-r_imaging*scale,'-.g','LineWidth',1.5);

%Exit pupil
p11=plot(z_e*ones(1,Ndraw)*scale,rho_exup*scale,'b','LineWidth',2);
plot(z_e*ones(1,Ndraw)*scale,rho_exdown*scale,'b','LineWidth',2);

%Vertices
p3=plot(V1*ones(1,Ndraw)*scale,rsample*scale,'--r','LineWidth',1.5);
p4=plot(V2*ones(1,Ndraw)*scale,rsample*scale,'--','color','#b22929','LineWidth',1.5);
p5=plot(V3*ones(1,Ndraw)*scale,rsample*scale,'--','color','#0095c9','LineWidth',1.5);
p6=plot(d_eye*ones(1,Ndraw)*scale,rsample*scale,'--k','LineWidth',1.5);
p7=plot(20e-3*ones(1,Ndraw)*scale,rsample*scale,'--g','LineWidth',1.5);

%Contours
plot(z_eye*scale,rho_eye*scale,'k','LineWidth',2);
p14=plot(z_cornea*scale,rho_cornea*scale,'color','#a0dfd4','LineWidth',2);
p15=plot(z_aIOL*scale,rho_IOL*scale,'color','#b22929','LineWidth',2);
p16=plot(z_pIOL*scale,rho_IOL*scale,'color','#0095c9','LineWidth',2);

hold off
xlabel('\fontsize{25} z [mm]');ylabel('\fontsize{25} r [mm]');
set(gca,'FontSize',25)
legend([p1 p2 p3 p4 p5 p6 p7 p11 p14 p15 p16],{'Retinal cone','Imaging cone','V_{c}','V_{a}','V_{p}','z_{r}','z_{i}','EP','LC','AI','PI'},'location','best');
grid on
xlim([0 25]);

figure(3)
%Close up on lens
hold on

%Rays
% p1=plot(z_rays*scale,r_rays*scale,'-.k','LineWidth',1.5);
% p2=plot(z_rays*scale,-r_rays*scale,'-.k','LineWidth',1.5);

%Vertices
p3=plot(V1*ones(1,Ndraw)*scale,rsample*scale,'--r','LineWidth',1.5);
p4=plot(V2*ones(1,Ndraw)*scale,rsample*scale,'--','color','#b22929','LineWidth',1.5);
p5=plot(V3*ones(1,Ndraw)*scale,rsample*scale,'--','color','#0095c9','LineWidth',1.5);
p6=plot(d_eye*ones(1,Ndraw)*scale,rsample*scale,'--k','LineWidth',1.5);
p17=plot(z_e*ones(1,Ndraw)*scale,rsample*scale,'r','LineWidth',2);


%Principal planes
p7=plot(H1*ones(1,Ndraw)*scale,rsample*scale,'m','LineWidth',1.5);
p8=plot(H2*ones(1,Ndraw)*scale,rsample*scale,'color','#a0dfd4','LineWidth',1.5);

%Aperture & Exit pupil
p9=plot(V2*ones(1,Ndraw)*scale,rho_apup*scale,'color','#FFA500','LineWidth',2);
p10=plot(V2*ones(1,Ndraw)*scale,rho_apdown*scale,'color','#FFA500','LineWidth',2);
p11=plot(z_e*ones(1,Ndraw)*scale,rho_exup*scale,'b','LineWidth',2);
p12=plot(z_e*ones(1,Ndraw)*scale,rho_exdown*scale,'b','LineWidth',2);

%Contours
p13=plot(z_eye*scale,rho_eye*scale,'k','LineWidth',2);
p14=plot(z_cornea*scale,rho_cornea*scale,'r','LineWidth',2);
p15=plot(z_aIOL*scale,rho_IOL*scale,'color','#b22929','LineWidth',2);
p16=plot(z_pIOL*scale,rho_IOL*scale,'color','#0095c9','LineWidth',2);

hold off
xlabel('\fontsize{25} z [mm]');ylabel('\fontsize{25} r [mm]');
legend([p4 p5 p17 p7 p8 p9 p11 p15 p16],{'V_{a}','V_{p}','DP','H_{1}','H_{2}','AS','EP','AI','PI'},'location','best');
set(gca,'FontSize',25)
grid on
xlim([3.5 4.0]);
ylim([-3.5 3.5]);