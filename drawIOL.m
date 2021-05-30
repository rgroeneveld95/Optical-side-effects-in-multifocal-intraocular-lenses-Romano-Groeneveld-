%Run master script
Optical_side_effects_in_multifocal_intraocular_lenses_mep; 
close all
clc



%Samples
Ndraw=100;
asample=linspace(0,2*pi,Ndraw);
scale=1000;

%Shape
t_bc=@(r) (abs(Rbc)-sqrt(Rbc^2-r.^2));                  %Base curve
t_ec=@(r,Rrad,Zrad) (Zrad-sqrt(Rrad^2-r.^2));           %Echelette curve
t_cen= @(r) abs(Rcen)-sqrt(Rcen^2-r.^2);                %Curvature of the central zone


%Interpolation of points between boundaries - Thickness functions
rint=zeros(length(r_i)-1,Ndraw);
thickbc=rint;
thickec=rint;
dthick=rint;
for i=1:length(r_i)-1
rint(i,:)=linspace(r_i(i),r_i(i+1),Ndraw);
thickbc(i,:)=t_bc(rint(i,:));
thickec(i,:)=t_ec(rint(i,:),Rrad(i),Zrad(i));
dthick(i,:)=thickec(i,:)-thickbc(i,:);
end
rint=rint'; rint=rint(:)';
thickbc=thickbc'; thickbc=thickbc(:)';
thickec=thickec'; thickec=thickec(:)';
dthick=dthick'; dthick=dthick(:)';

Rint=[linspace(0,r_i(1),Ndraw) rint linspace(r_i(end),r_iol,Ndraw)];
Thickbc2=t_bc(Rint);
Thickbc=[t_bc(linspace(0,r_i(1),Ndraw))+t_cen(linspace(0,r_i(1),Ndraw)) thickbc t_bc(linspace(r_i(end),r_iol,Ndraw))];
Thickec=[t_bc(linspace(0,r_i(1),Ndraw)) thickec t_bc(linspace(r_i(end),r_iol,Ndraw))];


%Visualization of the IOL lens

%2D
rlens=[fliplr(-Rint) Rint];
zant=[fliplr(Thickec) Thickec];
zpos=[fliplr(-Thickbc2) -Thickbc2]+2*Thickbc2(end);
r1=[fliplr(-linspace(0,r_i(1),Ndraw)) linspace(0,r_i(1),Ndraw)];
rlneg=fliplr(-linspace(r_i(end),r_iol,Ndraw));
rlpos=linspace(r_i(end),r_iol,Ndraw);
z1=t_bc(r1);    zlneg=t_bc(rlneg);  zlpos=t_bc(rlpos);

%3D
zant3d=repmat(Thickec,length(Thickec),1)';
zpos3d=repmat(-Thickbc+2*Thickbc(end),length(Thickbc),1)';
theta3d=linspace(0,2*pi,length(Rint));
x3d=zeros(length(Rint),length(Rint));
y3d=x3d;
for w=1:length(Rint)
x3d(w,:)=Rint(w)*cos(theta3d);
y3d(w,:)=Rint(w)*sin(theta3d);
end


%Plotting

% Front view IOL
figure(1)
hold on
for ap=2:i
    plot(r_i(ap)*cos(asample)*scale,r_i(ap)*sin(asample)*scale,'r','Linewidth',2);
end
plot(r_i(1)*cos(asample)*scale,r_i(1)*sin(asample)*scale,'b','Linewidth',2)
plot(r_iol*cos(asample)*scale,r_iol*sin(asample)*scale,'b','Linewidth',2)
hold off
grid on
title('\fontsize{25} Anterior surface view');xlabel('\fontsize{25} x [mm]');ylabel('\fontsize{25} y [mm]');
set(gca,'FontSize',25)
axis square


%Cross-section 2D view of the IOL
figure(2)
hold on
plot(zant*scale,rlens*scale,'r','Linewidth',2);
plot(zpos*scale,rlens*scale,'b','Linewidth',2);
plot(z1*scale,r1*scale,'b','Linewidth',2);
plot(zlneg*scale,rlneg*scale,'b','Linewidth',2);
plot(zlpos*scale,rlpos*scale,'b','Linewidth',2);
grid on
title('\fontsize{25} Cross-section');xlabel('\fontsize{25} z [mm]');ylabel('\fontsize{25} r [mm]');
set(gca,'FontSize',25)
axis tight

%3D view of the IOL
figure(3)
hold on
surf(zant3d*scale,x3d*scale,y3d*scale);
surf(zpos3d*scale,x3d*scale,y3d*scale);
camlight(-80,-45)
lighting phong
hold off
view(-45,10)
colormap('bone')
shading interp
view(-40,10)
camlight(-30,-45)
lighting phong
grid on
title('\fontsize{25} 3D view');xlabel('\fontsize{25} z [mm]');ylabel('\fontsize{25} y [mm]');zlabel('\fontsize{25} x [mm]');
set(gca,'FontSize',25)
axis square



% Overview of lens thickness

figure(4)
subplot(2,1,1)
hold on
ptec=plot(thickec*scale^2,rint*scale,'r','Linewidth',2);
ptbc=plot(thickbc*scale^2,rint*scale,'--b','Linewidth',2);
hold off
grid on
xlim([min(thickec*scale^2) max(thickec*scale^2)]);
xlabel('\fontsize{25} z [\mum]'); ylabel('\fontsize{25} r [mm]');
title('\fontsize{25} Surface profile diffractive region');
legend([ptbc ptec],{'\fontsize{25} Base curve','\fontsize{25} Echelette curve'},'location','best');
set(gca,'FontSize',25)


subplot(2,2,3)
hold on
plot(S_i*scale^2,r_i*scale,'--.','LineWidth',2,'MarkerSize',30)
xlabel('\fontsize{25} z [\mum]'); ylabel('\fontsize{25} r [mm]');
title('\fontsize{25} Diffractive step height profile');
set(gca,'FontSize',25)
grid on
hold off
axis tight

subplot(2,2,4)
hold on
plot((Thickec-Thickbc)*scale^2,Rint*scale,'Linewidth',2);
xlim([min((Thickec-Thickbc)*scale^2) max((Thickec-Thickbc)*scale^2)]);
xlabel('\fontsize{25} z [\mum]'); ylabel('\fontsize{25} r [mm]');
title('\fontsize{25} Relief profile optical region');
grid on
hold off
set(gca,'FontSize',25)


