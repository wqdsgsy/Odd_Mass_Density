clear all
close all
clc
%%%rho_y = 16277;

lambda = 37e9;
mu = 27.4e9;

% lambda = 2.2188e9;
% mu = 1.415e9;

rho_x = 16277;
rho_y = 16277;

% rat = 0;
% rho_xy = (rho_x+rho_y)/2*rat;

omega = 16e3*2*pi;
%
pp = 0;
%for omega1 = 5e3*2*pi:3e3*2*pi:20e3*2*pi
rat = 0.7;
     %rho_xy = (rho_x+rho_y)/2*rat;


%  pp = 0;
 for rat = 0.002:0.002:1
%  for theta = 0:pi/100:2*pi
rho_xy = rat*rho_x;
theta = pi/4;
    pp = pp + 1;
    M_l = [(lambda+2*mu)*(cos(theta))^2 + mu*(sin(theta))^2 .../
            (lambda+mu)*(cos(theta))*(sin(theta)); .../
            (lambda+mu)*(cos(theta))*(sin(theta)) .../
            (lambda+2*mu)*(sin(theta))^2 + mu*(cos(theta))^2];
       % omega = omega1;
    M_r = omega^2*[rho_x rho_xy; 0 rho_y];
    [Vec, Deig] = eig(M_r, M_l);
    theta_p(pp) = theta;
    rat_p(pp) = rat;
    if real(Deig(1,1)) < real(Deig(2,2))
        k2_r1(pp) = real(Deig(1,1));
        k2_r2(pp) = real(Deig(2,2));
        k2_i1(pp) = imag(Deig(1,1));
        k2_i2(pp) = imag(Deig(2,2));
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
        u_r1(pp) = real(Vec(1,1));
        v_r1(pp) = real(Vec(2,1));
        u_r2(pp) = real(Vec(2,1));
        v_r2(pp) = real(Vec(2,2));
        u_i1(pp) = imag(Vec(1,1));
        v_i1(pp) = imag(Vec(2,1));
        u_i2(pp) = imag(Vec(2,1));
        v_i2(pp) = imag(Vec(2,2));
        r1(pp) = real(Vec(2,1)/Vec(1,1));
        r2(pp) = real(Vec(2,2)/Vec(1,2));
        i1(pp) = imag(Vec(2,1)/Vec(1,1));
        i2(pp) = imag(Vec(2,2)/Vec(1,2));
        
        V1 = [Vec(1,1)/sqrt(abs(Vec(1,1))^2+abs(Vec(2,1))^2),Vec(2,1)/sqrt(abs(Vec(1,1))^2+abs(Vec(2,1))^2)];
        V2 = [Vec(1,2)/sqrt(abs(Vec(1,2))^2+abs(Vec(2,2))^2);Vec(2,2)/sqrt(abs(Vec(1,2))^2+abs(Vec(2,2))^2)];
        product(pp) = conj(V1)*V2;
    else
        k2_r2(pp) = real(Deig(1,1));
        k2_r1(pp) = real(Deig(2,2));
        k2_i2(pp) = imag(Deig(1,1));
        k2_i1(pp) = imag(Deig(2,2));
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
        u_r2(pp) = real(Vec(1,1));
        v_r2(pp) = real(Vec(2,1));
        u_r1(pp) = real(Vec(2,1));
        v_r1(pp) = real(Vec(2,2));
        u_i2(pp) = imag(Vec(1,1));
        v_i2(pp) = imag(Vec(2,1));
        u_i1(pp) = imag(Vec(2,1));
        v_i1(pp) = imag(Vec(2,2));
        r1(pp) = real(Vec(2,1)/Vec(1,1));
        r2(pp) = real(Vec(2,2)/Vec(1,2));
        i1(pp) = imag(Vec(2,1)/Vec(1,1));
        i2(pp) = imag(Vec(2,2)/Vec(1,2));
        V1 = [Vec(1,1)/sqrt(abs(Vec(1,1))^2+abs(Vec(2,1))^2),Vec(2,1)/sqrt(abs(Vec(1,1))^2+abs(Vec(2,1))^2)];
        V2 = [Vec(1,2)/sqrt(abs(Vec(1,2))^2+abs(Vec(2,2))^2);Vec(2,2)/sqrt(abs(Vec(1,2))^2+abs(Vec(2,2))^2)];
        product(pp) = conj(V1)*V2;
    end
    
end

figure(1);
% 

% end
% figure(1);
% % 
% plot(k_x_r1, k_y_r1)
% hold on
% plot(k_x_r2, k_y_r2)
% hold on;
% axis equal 
% xlabel('Re(k_x)');
% ylabel('Re(k_y)');

% plot(theta_p, k_i1)
% hold on
% plot(theta_p, k_i2)

plot(rat_p, k_r1,'r',rat_p, k_r2,'b:','linewidth',4);
set(gca, 'fontsize', 20);
xlabel('$\tilde{\beta}$','FontName', 'Times New Roman','Interpreter','latex', 'fontsize', 25);
ylabel('$\Re(q)$ (1/m)','FontName', 'Times New Roman','Interpreter','latex', 'fontsize', 25);

%grid on;
set(gca, 'FontName', 'Times New Roman');
set(gca, 'linewidth',1.5);
  figure(2);
% % 
%  plot(k_x_i1, k_y_i1)
% hold on
%  plot(k_x_i2, k_y_i2)
% axis equal 
% xlabel('Im(k_x)');
% ylabel('Im(k_y)');

% plot(theta_p, k_i1)
% hold on
% plot(theta_p, k_i2)
% xlabel('\theta');
% ylabel('Im(k)');

plot(rat_p,k_i1,'r',rat_p,k_i2,'b:','linewidth',4);
set(gca, 'fontsize', 20);
xlabel('$\tilde{\beta}$','FontName', 'Times New Roman','Interpreter','latex', 'fontsize', 25);
ylabel('$\Im(q)$ (1/m)','FontName', 'Times New Roman','Interpreter','latex', 'fontsize', 25);

set(gca, 'linewidth',1.5);
%grid on;
set(gca, 'FontName', 'Times New Roman');
 ratt = 0.002:0.002:1;
figure(3);
plot(ratt,r1,'r',ratt,r2,'b:','linewidth',4);
ylim([-10,10]);
set(gca, 'fontsize', 20);
xlabel('$\tilde{\beta}$','FontName', 'Times New Roman','Interpreter','latex', 'fontsize', 25);
ylabel('$\Re(\tilde{U}_2/\tilde{U}_1)$','FontName', 'Times New Roman','Interpreter','latex', 'fontsize', 25);

set(gca, 'FontName', 'Times New Roman');
set(gca, 'linewidth',1.5);
%grid on;
%set(gca,'XMinorTick','on','YMinorTick','on');

figure(4);
plot(ratt,i1,'r',ratt,i2,'b:','linewidth',4);
set(gca, 'fontsize', 20);
xlabel('$\tilde{\beta}$','FontName', 'Times New Roman','Interpreter','latex', 'fontsize', 25);
ylabel('$\Im(\tilde{U}_2/\tilde{U}_1)$','FontName', 'Times New Roman','Interpreter','latex', 'fontsize', 25);

set(gca, 'FontName', 'Times New Roman');
set(gca, 'linewidth',1.5);

%grid on;
%set(gca,'XMinorTick','on','YMinorTick','on');
% figure(5);
% plot(ratt,product);
% 
