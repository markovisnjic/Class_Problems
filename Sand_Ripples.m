%%Sand Ripple Evolution 1D
%Marko Visnjic

clf
figure(1)

%%Parameters

beta = .1;
S = 10^7; %m^-2/s^-1

%% Arrays

D = 0.005; % grain diameter
dx = 5; 
xmax = 1000;
x =0:dx:xmax;

g_remove = 10; %Grains removed with impact
g_hop = 5*dx; %How far grains go after ejection

z0=0.1;
%z = z0+0.01*rand(size(x)); % random topography with a st dev 0.01
z = z0*ones(size(x));
porosity = 0.35;
N = ceil(dx*z*(1-porosity)/(pi*D*D/4)); % number of grains in each bin.
z = (pi*N*D^2)/(4*(1-porosity)*dx);

% Number of Impacts
dt = 1;
tmax = 100000;
t = 0:dt:tmax;
imax =length(t);


%Particle Impacts
for i=1:imax
    
particle_wind = 0.01*((xmax)*rand-x); %generate the particle trajectories

intercept = find(particle_wind < z); %find locations where the particle wind impacts the topo
  
Impact_Location = intercept(1); % take the first index as the impact point

   
if (intercept > (length(N)-g_hop)) %Wrap around B.C.
    N(Impact_Location) = N(Impact_Location)-g_remove; %Removing on impact
    g_enter = g_hop-(length(N)-Impact_Location);
    N(g_enter) = N(g_enter)+g_remove;
else
     N(Impact_Location) = N(Impact_Location)-g_remove; %Removing on impact
     N(Impact_Location+g_hop) = N(Impact_Location+g_hop)+g_remove; %Add grains

end
z = (pi*N*D^2)/(4*(1-porosity)*dx);
z = max(z0,z);
hold on
%plot(x,particle_wind)
plot(x/100,z,'k')
hold off
ylabel('Height (m)','FontSize',24,'fontname','arial')
xlabel('Distance (m)','FontSize',24,'fontname','arial')

%pause(0.1)
drawnow
%end
end