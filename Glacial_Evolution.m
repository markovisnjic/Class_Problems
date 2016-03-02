%% Glacial Evolution
% Marko Visnjic 2/12/2016
% (Worked with Kelly Kochanski)
%% Parameters
clf
figure (1)
rho_i = 917; %Density of Ice
gamma = 0.01; %y^-1
s = 0.4;
g = 9.81;
A = 2.1*10^-16;
%% Arrays
dx = 1;
xmax = 1000;
x = 0:dx:xmax; %Distance Array m

dz = 0.1;
zmax = 100;
zb =zmax:-dz:0; %Elevation Array m
ELA = (2/3)*zmax;

b = gamma*(zb-ELA);
G0 = 30+30.*(-(x-(xmax/3)).^2)/(xmax/3).^2; %Initial Glacial thickness m
G0 = max(0, G0);
G = G0.*ones(size(x)); %Initial Glacial Array

%Elevation of glacier surface
h = (zmax-s*x)+G; %m

%Time
dt = 1/52;
tmax = 1000;
t = 0:dt:tmax;
imax = length(t);

nplots = 100;
tplot = tmax/nplots;


for i =1:imax
    r_num(i) = rand(1);
   % array of fluxes of ice
   q = A*(rho_i*g*s).^3.*G.^5/5;
   q = [0 q];
   dqdx = diff(q)/dx;
   ELA = (2/3)*zmax+20*r_num(i);
   dGdt = b - dqdx;
   % update glacier thickness and elevation
   G = G+dGdt*dt;
   G = max(G, 0);
   h = G + zb;
   % plot all the things every 100 steps
   if rem(t(i),tplot)==0
    plot(x,h,'b')
    hold on
    plot(x,zb,'r--')
    plot(x,ELA*ones(size(x)),'k--')
    axis([0 x(end) 0 150])
      xlabel('Distance','fontname','arial','fontsize',21)
        ylabel('Elevation','fontname','arial','fontsize',21)
        set(gca,'fontsize',18,'fontname','arial')
         title ('Glacial Evolution')
    drawnow
   end
end