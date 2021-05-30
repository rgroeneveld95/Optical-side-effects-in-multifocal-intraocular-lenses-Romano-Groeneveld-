%Run master script
Optical_side_effects_in_multifocal_intraocular_lenses_mep; 
close all
clc


%DATA
rpupdata=[1 2 3 4 5 6];
If3=[38.54 44.711 57.77 74.93 84.00 88.87];
In3=[42.51 37.92 28.74 16.40 9.57 7.21];

If1=[38.78 46.01 58.92 75.64 84.44 89.17];
In1=[43.23 37.89 28.95 16.48 10.46 7.23];



%Samples
scale=1000;
Ndraw=1000;
r_pupil=3e-3;

order=0:4;
a_i=(n_iol-n_a)/l_d*S_i;   %Induced phase shifts (in units of wavelength)
area=zeros(1,length(r_i)-1);
area(1)=pi*r_i(1)^2;
for surface=2:length(r_i)
area(surface)=pi*(r_i(surface)^2-r_i(surface-1)^2);
end
t_e=zeros(length(order),length(r_i));
for m=1:length(order)
    m_th_order=order(m);
t_e(m,:)=sinc(m_th_order-a_i).^2;
end
energy=t_e.*area;
energytot=sum(energy,2);

boundary=[0 r_i r_pupil];
t_epupil=[t_e [1; zeros(1,length(order)-1)']];
area_pupil=[area pi*(r_pupil.^2-r_i(end).^2)];
e_pupil=t_epupil.*area_pupil;
pr=linspace(0.5e-3,r_pupil,Ndraw);


e_tot=sum(sum(e_pupil));
E_pupil=cell(length(order),length(boundary)-1);
for m=1:length(order)
    for z=2:length(boundary)
        if z==2
        E_pupil{m,z-1}=@(r) ((e_pupil(m,z-1).*(r.^2-boundary(z-1).^2)*pi./area_pupil(z-1)).*(r>boundary(z-1) & r<=boundary(z)))./(pi*r.^2);
        else
        E_pupil{m,z-1}=@(r) ((sum(e_pupil(m,1:z-2))+e_pupil(m,z-1).*(r.^2-boundary(z-1).^2)*pi./area_pupil(z-1)).*(r>boundary(z-1) & r<=boundary(z)))./(pi*r.^2);
        end
    end
end

Epupil=zeros(length(order),Ndraw);
for m=1:length(order)
for z=1:length(boundary)-1
Epupil(m,:)=Epupil(m,:)+sum(E_pupil{m,z}(pr),1);
end 
end



%Plotting

%Light distribution

figure(1)
hold on
p1=plot(pr*2*scale,Epupil(1,:)*100,'--r','Linewidth',2);
p2=plot(pr*2*scale,Epupil(2,:)*100,'--b','Linewidth',2);
% p3=plot(rpupdata,If1,'.-r','LineWidth',2,'MarkerSize',25);
% p4=plot(rpupdata,In1,'.-b','LineWidth',2,'MarkerSize',25);
p3=plot(rpupdata,If3,'.-r','LineWidth',2,'MarkerSize',25);
p4=plot(rpupdata,In3,'.-b','LineWidth',2,'MarkerSize',25);

hold off
xlim([min(2*pr*scale) max(2*pr*scale)+1e-20]);
legend([p1 p2 p3 p4 ],{'\fontsize{25} Theoretical far focus','\fontsize{25} Theoretical near focus','\fontsize{25} Simulated near focus','\fontsize{25} Simulated far focus'},'Location','Best')
% legend('\fontsize{25} Distant','\fontsize{25} Near','Location','East');
xlabel('\fontsize{25} IOL Pupil diameter [mm]');
ylabel('\fontsize{25} Energy efficiiency [%]');
set(gca,'FontSize',25)
grid on
% axis square