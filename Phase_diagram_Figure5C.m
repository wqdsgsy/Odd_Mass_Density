clear all
close all
clc

lambda = 37e9;
mu = 27.4e9;

E = 70.54e9;
nu = 0.2873;

gamma = lambda/mu;

theta1 = linspace(0,2*pi,2001);
beta1 = linspace(-5,5,2001);
for i = 1:2001
theta = theta1(i); 
for j = 1:2001
    m = sin(theta);
    n = cos(theta);
   beta = beta1(j);
   A = E/(1-nu^2)*cos(theta)^2+E/2/(1+nu)*sin(theta)^2;
   C = E/(1-nu^2)*sin(theta)^2+E/2/(1+nu)*cos(theta)^2;
   B = E/2/(1-nu)*cos(theta)*sin(theta);
   ref=beta^2*B^2-2*B*(A+C)*beta+4*B^2+(A-C)^2;
   if ref<0
    AA(j,i) = 0;
   else if beta>6 && theta>20/180*pi && theta<70/180*pi 
    AA(j,i) = 0.5;
       else if beta>6 && theta>190/180*pi && theta<260/180*pi 
    AA(j,i) = 0.5;
           else if beta<-6 && theta>280/180*pi && theta<350/180*pi 
    AA(j,i) = 0.5;
       else if beta<-6 && theta>110/180*pi && theta<160/180*pi 
    AA(j,i) = 0.5;
       else
           AA(j,i) = 1;
   end
               end
           end
       end
end
end
end
[t,r] = meshgrid(linspace(0,2*pi,2001),linspace(-5,5,2001));

[theta1,beta1] = pol2cart(t,r); 
%[R PHI] = meshgrid(beta1 ,theta1);
%h1=surf((R+5).*cos(PHI), (R+5).*sin(PHI), AA);
figure('color','white');  
polarplot3d(AA,'PlotType','surfn','PolarGrid',{4 8},'TickSpacing',22.5,...
    'AngularRange',[0 2]*pi,'RadialRange',[0 5]);

set(gca,'DataAspectRatio',[1 1 10],'View',[-12,38],...
    'Xlim',[-5.5 5.5],'Xtick',[-5 0 5],...
    'Ylim',[-5.5 5.5],'Ytick',[-5 0 5]);

view(0,90);
set(gca,'XMinorTick','off','YMinorTick','off');
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
set(gca,'fontsize',12);
set(gca,'XColor', 'none','YColor','none');

print(gcf,'-r1200','-dpng','./phase dia');