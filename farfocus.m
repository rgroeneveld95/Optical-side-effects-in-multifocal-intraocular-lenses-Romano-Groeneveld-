clear
close all
clc

%=============================================================================================
%OPTICAL SIDE EFFECTS IN MULTIFOCAL INTRAOCULAR LENSES
%=============================================================================================

%SYMBOL EXPLANATION

%Vertices

%V1 = the vertex of the cornea
%V2 = the vertex of the anterior surface of the smooth IOL lens
%V3 = the vertex of the posterior surface of the smooth IOL lens

%Characteristic lengths

%so1 = object distance w.r.t. primary principal plane H1
%T1 = distance of the primary principal plane H1 to the primary vertex V1
%d = anterior chamber width measured from vertex V1 to vertex V2
%d_iol = width of the IOL measured from vertex V2 to vertex V3
%T2 = distance of the secondary principal plane H2 to the last vertex V3
%si2 = image distance w.r.t. secondary principal plane H2
%d_eye = the diameter from the eye
%r_c = the maximal radial height of the cornea

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
%LOCATING INPUT PARAMETERS
%=============================================================================================

%Wave field parameters
lambda=550e-9;                          %Wavelength of optical light-green color 
k=2*pi/lambda;                          %Wave number

%IOL lens input parameters
%(ReSTOR Aspheric IOL SN6AD3)
l_d=550e-9;                             %Design wavelength
Rbc=28.22e-3;                           %Base-curve radius of curvature
i=12;                                   %Amount of zones
i_in=1;                                 %Inner zone of apodization
i_out=12;                               %Outer zone of apodization
D_add=4;                                %Add power in dpt
p=0.5;                                  %Phase height [between 1.0 - 0.0]
ap_constant=3;                          %Apodization constant
r_iol=3.0e-3;                           %Maximal radial height of the IOL
r_cen=0.742e-3/2;                       %Maximal radial height of the central zone of the IOL
d_in=0e-3;                              %The inner thickness of the IOl, excluding radial surfaces

%IOL lens input parameters
%(ReSTOR Aspheric IOL SN6AD1)
% l_d=550e-9;                             %Design wavelength
% Rbc=28.22e-3;                           %Base-curve radius of curvature
% i=9;                                    %Amount of zones
% i_in=1;                                 %Inner zone of apodization
% i_out=9;                                %Outer zone of apodization
% D_add=3;                                %Add power in dpt
% p=0.5;                                  %Phase height [between 1.0 - 0.0]
% ap_constant=3;                          %Apodization constant
% r_iol=3.0e-3;                           %Maximal radial height of the IOL
% r_cen=0.856e-3/2;                       %Maximal radial height of the central zone of the IOL
% d_in=0e-3;                              %The inner thickness of the IOl, excluding radial surfaces

%Refractive indices
n_iol=1.5542;
n_air=1.00;                             
n_a=1.336;                              
n_v=1.336;

%Radii of curvature
Rc=7.645e-3;                              
Ra=Rbc;                                
Rp=-Rbc;                    

%Characteristic lengths
d=3.6e-3;
d_iol=2*(Rbc-sqrt(Rbc.^2-r_iol.^2))+d_in;
d_eye=24.17e-3;
r_c=11.5e-3/2;
r_pupil=3.0e-3;
s01=40;                 %Far focus distance
% s01=-0.3258;          %Near defocus distance SN6AD3
% s01=-0.4350;          %Near defocus distance SN6AD1
    
%Vertices
V1=0;
V2=V1+d;
V3=V2+d_iol;

%=============================================================================================
%GEOMETRICAL OPTICS CALCULATIONS
%=============================================================================================

%Calculate image and principal planes
[T1,T2,si2]=geo_image(s01,n_air,n_a,n_iol,n_v,Rc,Ra,Rp,d_iol,V1,V2);
H1=V2+T1;
H2=V3+T2;

%Calculate exit pupil location
s_aperture=T1;
s_exit=exit_pupil(s_aperture,n_a,n_iol,n_v,Ra,Rp,d_iol);

%Magnifications
Mt=abs(s_exit/s_aperture);
Ml=Mt^2;

%Location of geometrical image
z_im=H2+si2;

%Location of the aperture
z_a=V2;

%Exit pupil location and radial height
z_e=H2+s_exit;
r_e=Mt*r_pupil;

%Distance geometric image w.r.t. the exit pupil
d_i=z_im-z_e;

%Distance retinal image plane w.r.t the exit pupil
d_r=d_eye-z_e;

%=============================================================================================
%IOL LENS STRUCTURE
%=============================================================================================

%Radial boundaries of the anterior surface of the IOL
r_i=sqrt((2*(1:i)-1)*l_d/D_add);        %Boundary points of the annular zones
r_i(1)=r_cen;                           %Set 1st annular boundary to special radial length
r_in=r_i(i_in);                         %Starting boundary point of apodization zone
r_out=r_i(i_out);                       %Ending boundary point of apodization zone


% General Stepheight and Apodization

%Stepheight
% S=l_d/(2*(n_iol-n_a));
% %Apodization
% f_ap=zeros(1,length(r_i));                                                        
% for ap=1:i
%     if r_i(ap)<r_in                                                   %If r_i <= r_in
%         f_ap(ap)=1;
%     elseif r_i(ap)>r_out                                              %If r_i >= r_out
%         f_ap(ap)=0;
%     else
%         f_ap(ap)=1-((r_i(ap)-r_in)/(r_out-r_in)).^ap_constant;        %If r_in <r_i <= r_out
%     end
% end


% Stepheight and Apodization of reSTOR SN6AD3

%Stepheight
S=1.31189e-6*1; 
%Apodization
f_ap=1-(r_i/r_i(end)).^ap_constant;


%Reduced stepheight
S_i=S*f_ap;      

%Z-coordinate boundaries of the anterior surface of the IOL

%z-coordinate of the boundary points on basecurve
z_i=abs(Rbc)-sqrt(Rbc^2-r_i.^2);

%z_coordinate of the boundary points on echelette curve
z_e_down=z_i-p*S_i;             %stepheight beneath base curve
z_e_up=z_i+(1-p)*S_i;           %stepheight above base curve

%Sphere paramaters of the echelette curve (the centres and radii of curvature)
zone=i-1;                       %No. surfaces= No.zones -1
Zrad=zeros(1,zone);
Rrad=Zrad;
for i=1:zone
Zrad(i)=(r_i(i+1)^2-r_i(i)^2+z_e_up(i+1)^2-z_e_down(i)^2)/(2*(z_e_up(i+1)-z_e_down(i)));
Rrad(i)=sqrt((z_e_up(i+1)-Zrad(i))^2+r_i(i+1)^2);
end

%Radii of curvature of the central zone
% correction=+20e-3;
correction=0;
Rcen=(r_i(1)^2+(S_i(1)/2)^2)/(2*S_i(1)/2)+correction;  

%Radial boundaries
r_b=[0 r_i r_iol];

%Thickness functions
t_bc=@(r) (abs(Rbc)-sqrt(Rbc^2-r.^2));                  %Base curve
t_ec=@(r,Rrad,Zrad) (Zrad-sqrt(Rrad^2-r.^2));           %Echelette curve
t_cen= @(r) abs(Rcen)-sqrt(Rcen^2-r.^2);                %Curvature of the central zone
% t_cen= @(r) 0;                %Curvature of the central zone



%Depth structure of the diffraction plate
dz=cell(1,length(r_b)-3);
for b=3:length(r_b)-1
dz{b-2}=@(r) (t_ec(r,Rrad(b-2),Zrad(b-2))-t_bc(r)).*(r>r_b(b-1) & r<r_b(b));
end
DZ=@(r) (sum(cell2mat(cellfun(@(fun)fun(r),dz(:),'uni',0)),1)+(t_cen((r))).*(r>0 & r<r_b(2)));
DZ=@(r,C) DZ(r)*C;

%=============================================================================================
%DIFFRACTION INTEGRAL FOR APERTURE AND DIFFRACTION PLATE
%=============================================================================================

%Induced phase due to the diffraction plate including appropiate magnifications
T=@(r,C) exp(1i*k*(n_iol-n_a).*DZ(r/Mt,C)*Ml)/Mt;

%Prefactor in front of the fresnel integral
prefac=@(rv) 1i*2*pi*n_v*n_air./(lambda^2*s01*d_r).*exp(1i*k*n_v./(2*d_r).*rv.^2);

%Quadtratic phase term in the in fresnel integral
phase=@(r,d_i) exp(1i*k*n_v/2.*(1/d_r-1./d_i).*r.^2);

%Bessel function of order zero
J0=@(r,rv) besselj(0,k*n_v./d_r.*r.*rv);

%Pupil function of the exit pupil
pupil=@(r) (r>0 & r<r_e);

%Integrand of the diffraction integral
I=@(r,rv,d_i,C) pupil(r).*J0(r,rv).*phase(r,d_i).*T(r,C).*r;


%Sampling exit pupil plane coordinates and image plane coordinates
N=1000;
r=linspace(0,3e-3*Mt,N);
rv=linspace(0,15e-6,N);

%Aperture
U_i0=zeros(1,N);
for n=1:N
U_i0(n)=prefac(rv(n)).*trapz(r,I(r,rv(n),d_i,0));
end
I_i0=abs(U_i0).^2;
% I_i0=abs(U_i0).^2*s01^2;
P_i0=2*pi*trapz(rv,I_i0.*rv);

%Diffraction plate
U_i1=zeros(1,N);
for n=1:N
U_i1(n)=prefac(rv(n)).*trapz(r,I(r,rv(n),d_i,1));
end
I_i1=abs(U_i1).^2;
% I_i1=abs(U_i1).^2*s01^2;
P_i1=2*pi*trapz(rv,I_i1.*rv);

%=============================================================================================
%FINDING FOCAL POINTS
%=============================================================================================

d_im=linspace(16e-3,26e-3,N);
s_ob=geo_image_inverse(d_im+z_e-H2,n_air,n_a,n_iol,n_v,Rc,Ra,Rp,d_iol,V1,V2);

s_ob2=linspace(-30,30,N);
[T1,T2,s_im]=geo_image(s_ob2,n_air,n_a,n_iol,n_v,Rc,Ra,Rp,d_iol,V1,V2);
d_im2=s_im+H2-z_e;

%Aperture (Should be 1 clear focus)
%Calculate the field for rv=0 (only on the optical axis)
U_im0=zeros(1,N);
for n=1:N
    U_im0(n)=prefac(0).*trapz(r,I(r,0,d_im(n),0));
% U_im0(n)=prefac(0).*trapz(r,I(r,0,d_im(n),0))/r_e^2;
end
I_im0=abs(U_im0).^2;

%Diffraction plate (Should give more focii)
%Calculate the field for rv=0 (only on the optical axis)
U_im1=zeros(1,N);
for n=1:N
    U_im1(n)=prefac(0).*trapz(r,I(r,0,d_im(n),1));
%    U_im1(n)=prefac(0).*trapz(r,I(r,0,d_im(n),1))/r_e^2;
end
I_im1=abs(U_im1).^2;

%=============================================================================================
%VISUALIZATION AND PLOTS
%=============================================================================================

%Plotting results of the diffraction integrals
rv_plot=[-fliplr(rv) rv(2:end)];
I_i0_plot=[fliplr(I_i0) I_i0(2:end)];
I_i1_plot=[fliplr(I_i1) I_i1(2:end)];

%Analytic solution for aperture
PSF=(n_air*r_e./(lambda*s01)).*besselj(1,k*n_v*r_e.*rv/d_r)./rv;
airypattern=abs(PSF).^2;
airyplot=[fliplr(airypattern) airypattern(2:end)];


%PSF Results (1D)
figure(1)
hold on
plot(rv_plot*10^6,airyplot,'r','LineWidth',2)
semilogy(rv_plot*10^6,I_i0_plot,'--b','LineWidth',2)
% semilogy(rv_plot*10^6,I_i1_plot,'r','LineWidth',2)
xlabel('\fontsize{25} r [\mum]');ylabel('\fontsize{25} Intensity');
% legend({['\fontsize{25} Aperture P=',num2str(P_i0,4)],['\fontsize{25} Diffraction plate P=',num2str(P_i1,4)]},'Location','northeast')
legend('\fontsize{25} J_{1}','\fontsize{25} IOL','Location','NorthEast')
hold off
set(gca,'FontSize',25)
grid on
axis square

%PSF Results (2D)
figure(2)

I_i1_2d=repmat(I_i1,N,1)';
theta2d=linspace(0,2*pi,N);
x2d=zeros(N,N);
y2d=x2d;
for w=1:N
x2d(w,:)=rv(w)*cos(theta2d);
y2d(w,:)=rv(w)*sin(theta2d);
end

hold on
mesh(x2d*10^6,y2d*10^6,I_i1_2d);
view(0,90);
xlim([-max(rv)/sqrt(2) max(rv)/sqrt(2)]*10^6);
ylim([-max(rv)/sqrt(2) max(rv)/sqrt(2)]*10^6);
colorbar 
hold off
axis square
xlabel('\fontsize{25} x [\mum]');
ylabel('\fontsize{25} y [\mum]');
set(gca,'FontSize',25)


%FOCAL POINTS
figure(3)
hold on
% plot((d_r-d_im)*1000,I_im0,'b')
% plot((d_r-d_im)*1000,I_im1,'r')
plot((d_im)*1000,I_im0,'b')
plot((d_im)*1000,I_im1,'r')
% xlabel('z_i - z_r [mm]'); 
xlabel('z_i [mm]'); 
ylabel('Center intensity');
title(['Center intensity as function image distance for r_{pupil}=',num2str(r_pupil*1000,3),' mm'])
legend('Aperture','Diffraction plate')
hold off


%PROPAGATION-OBJECT DISTANCE RELATION
figure(4)
hold on
plot(s_ob2,d_im2*1000,'b','Linewidth',2)
xlim([min(s_ob2) max(s_ob2)]);
xlabel('\fontsize{25} S_{1o} [m]');
ylabel('\fontsize{25} z_{i} [mm]');
set(gca,'FontSize',25)
grid on
% plot(s_ob,d_r*ones(1,length(s_ob))*1000,'r')
hold off
