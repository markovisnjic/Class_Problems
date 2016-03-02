%% Glacial Evolution on Mountain Top
% Marko Visnjic 2/12/2016
%% Parameters
figure (1)
clf

rho_i = 917.0; %Density of Ice
gamma = 0.01; %y^-1
s = 0.05;
g = 9.81;
A = 2.1*10^-16;
%% Arrays

% Mountain

dxM = 200;
xMmax = 40000;
xM = -xMmax:dxM:xMmax;

%Elevation
dzM = 1;
zMmax = 4000;
zM = zMmax:-dzM:0;


%% Building the Glacier
ELA = 3600;
b = gamma*(zM-ELA);
G = zeros(size(xM));

%Initial Topography

zB = zMmax-(s*abs(xM))+(xMmax/500)*cos((xM)/0.35)-((xMmax/1000)*sin((xM)/3))-(xMmax/2000)*sin(3*(xM)/5);%Mountain top geography

h = zB+G; %Total height 

%%Time
dt = 0.002;
tmax = 1000;
t = 0:dt:tmax;
imax = length(t);
nplots = 40;
tplot=tmax/nplots;


for i = 1:imax
    b = gamma*(h-ELA);
    slope = diff(h)/dxM;
    Gedge = G(1:end-1)+(0.5*diff(G));
  q = -sign(slope).*(A*(rho_i*g*abs(slope)).^3).*(Gedge.^5)/5; %flux rule
  q = [0 q 0];
  dqdx = diff(q)/dxM;
  dGdt = b - dqdx;  
  G = G+(dGdt*dt);
  G = max(G,0);
  h = G + zB;
  
  if(rem(t(i),tplot)==0)
      figure(1)
  plot(xM/1000,h,'c')
  hold on
  plot(xM/1000,zB,'k','linewidth',2)
  set(gca,'fontsize',21,'fontname','arial');
ylabel('Elevation (m)','FontSize',24,'fontname','arial')
xlabel('Distance (km)','FontSize',24,'fontname','arial')

  pause(0.1)
  end
  
end