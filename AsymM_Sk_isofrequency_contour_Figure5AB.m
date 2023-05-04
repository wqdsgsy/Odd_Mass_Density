clear all
close all
clc
%
lambda = 37e9;
mu = 27.4e9;
%
rho_x = 16277;
rho_y = 16277;

rho_xy = 0.7*rho_x;
%rho_xy = 10000;
omega1 = linspace(5e3,20e3,16)*2*pi;
%
pp = 0;
% rho_y = 16277*4;
% for rat = 0:0.002:1
% rho_xy = (rho_x + rho_y)/2*rat;
% for i = 1:16
     omega = 12.5e3*2*pi;
for theta = 0:pi/100:2*pi
% for theta = pi/4
    pp = pp + 1;
    M_l = [(lambda+2*mu)*(cos(theta))^2 + mu*(sin(theta))^2 .../
            (lambda+mu)*(cos(theta))*(sin(theta)); .../
            (lambda+mu)*(cos(theta))*(sin(theta)) .../
            (lambda+2*mu)*(sin(theta))^2 + mu*(cos(theta))^2];
    M_r = omega^2*[rho_x rho_xy; 0 rho_y];
    M_try = inv(M_r)*M_l
    [Vec, Deig] = eig(M_r, M_l);
    theta_p(pp) = theta;
    %rat_p(pp) = rat;
    xx(pp) = Deig(1,1);
    yy(pp) = Deig(2,2);
     if real(Deig(1,1)) < real(Deig(2,2))
        k_x_r1(pp) = real(sqrt(Deig(1,1)))*cos(theta);
        k_y_r1(pp) = real(sqrt(Deig(1,1)))*sin(theta);
        k_x_r2(pp) = real(sqrt(Deig(2,2)))*cos(theta);
        k_y_r2(pp) = real(sqrt(Deig(2,2)))*sin(theta);
        k_x_i1(pp) = imag(sqrt(Deig(1,1)))*cos(theta);
        k_y_i1(pp) = imag(sqrt(Deig(1,1)))*sin(theta);
        k_x_i2(pp) = imag(sqrt(Deig(2,2)))*cos(theta);
        k_y_i2(pp) = imag(sqrt(Deig(2,2)))*sin(theta);
        k_r1(pp) = real(sqrt(Deig(1,1)));
        k_r2(pp) = real(sqrt(Deig(2,2)));
        k_i1(pp) = imag(sqrt(Deig(1,1)));
        k_i2(pp) = imag(sqrt(Deig(2,2)));
    else
        k_x_r2(pp) = real(sqrt(Deig(1,1)))*cos(theta);
        k_y_r2(pp) = real(sqrt(Deig(1,1)))*sin(theta);
        k_x_r1(pp) = real(sqrt(Deig(2,2)))*cos(theta);
        k_y_r1(pp) = real(sqrt(Deig(2,2)))*sin(theta);
        k_x_i2(pp) = imag(sqrt(Deig(1,1)))*cos(theta);
        k_y_i2(pp) = imag(sqrt(Deig(1,1)))*sin(theta);
        k_x_i1(pp) = imag(sqrt(Deig(2,2)))*cos(theta);
        k_y_i1(pp) = imag(sqrt(Deig(2,2)))*sin(theta);
        k_r2(pp) = real(sqrt(Deig(1,1)));
        k_r1(pp) = real(sqrt(Deig(2,2)));
        k_i2(pp) = imag(sqrt(Deig(1,1)));
        k_i1(pp) = imag(sqrt(Deig(2,2)));
    end
    
end
% end
figure(1);
plot(k_x_r1/max(k_y_r2),k_y_r1/max(k_y_r2),'r','linewidth',4)
hold on
plot(k_x_r2/max(k_y_r2), k_y_r2/max(k_y_r2),'b--','linewidth',4)
 xlim([-1,1]); ylim([-1,1]);
set(gca, 'linewidth',1.5);
pbaspect([1 1 1]);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
print(gcf,'-r600','-dpng','./ifc_a4');
%set(gca,'Yticklabel',[]) 
%set(gca,'Xticklabel',[]) 
%set(gca,'xtick',[])
%set(gca,'ytick',[])

% plot(k_x_i1, k_y_i1)
% hold on
% plot(k_x_i2, k_y_i2)

%%%%%%%%%%%%%%%%%
figure(2);
plot(theta_p, k_i1/max(k_i2),'r','linewidth',4)
hold on
plot(theta_p, k_i2/max(k_i2),'b--','linewidth',4)
pbaspect([1 1 1])
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
set(gca, 'linewidth',1.5);

xlim([0,2*pi]);
print(gcf,'-r600','-dpng','./ifc_b4');
%%%%%%%%%%%%%%%%%%

%ylim([-1,1]);
% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[]) 
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% plot(rat_p, k_i1)
% hold on
% plot(rat_p, k_i2)


