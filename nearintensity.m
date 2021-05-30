%PSF Results (1D)
figure(1)
hold on
% plot(rv_plot*10^6,airyplot,'g')
semilogy(rv_plot*10^6,I_i0_plot,'b','LineWidth',2)
semilogy(rv_plot*10^6,I_i1_plot,'r','LineWidth',2)
xlabel('\fontsize{25} r [\mum]');ylabel('\fontsize{25} Intensity');
% legend({['\fontsize{25} Aperture P=',num2str(P_i0,4)],['\fontsize{25} Diffraction plate P=',num2str(P_i1,4)]},'Location','northeast')
legend('\fontsize{25} J_{1}','\fontsize{25} SN6AD3','Location','NorthEast')
hold off
set(gca,'FontSize',25)
grid on
axis square
xlim([-15 15])

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
