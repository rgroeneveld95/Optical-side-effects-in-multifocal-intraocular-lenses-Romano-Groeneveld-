%Run master script
Optical_side_effects_in_multifocal_intraocular_lenses_mep; 
close all
clc

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

figure(1)
hold on
for pp=1:length(r_i)
plot([r_i(pp)*scale r_i(pp)*scale],[0 100],'--k','Linewidth',2)
txten={['\fontsize{25} r_{',num2str(pp)],'}'};
text(r_i(pp)*scale,98,txten,'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','color','black')
end
pp=plot(r_i*scale,t_e*100,'.-','Linewidth',2,'MarkerSize',30);
set(gca,'FontSize',25)
hold off
legend([pp(1) pp(2)],{'Distant','Near'},'Location','East');
grid on
xlabel('\fontsize{25} r [mm]');ylabel('\fontsize{25} Throughput efficiency [%]')
% axis square

%Light distribution

figure(2)
plot(pr*2*scale,Epupil*100,'Linewidth',2)
xlim([min(2*pr*scale) max(2*pr*scale)+1e-20]);
legend('\fontsize{25} Distant','\fontsize{25} Near','Location','East');
xlabel('\fontsize{25} IOL Pupil diameter [mm]');
ylabel('\fontsize{25} Energy efficiiency [%]');
set(gca,'FontSize',25)
grid on
% axis square