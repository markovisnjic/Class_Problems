%%Sand Ripple Evolution 1D
%Marko Visnjic

clf
figure(1)

%%Parameters

beta = .1;
S = 10^7; %m^-2/s^-1

%% Arrays

D = 0.005; % grain diameter
dx = 5; % make these bins much bigger than your grains
xmax = 1000;
x = 0:dx:xmax;

g_remove = 10; %Grains removed with impact
g_hop = 5*dx; %How far grains go after ejection

% here is how i propose you make your initial topo and 
% back out how many grains are in each bin
z0=1;
z = z0+0.01*rand(size(x)); % random topography with a st dev 0.01

porosity = 0.35;
N = ceil(dx*z*(1-porosity)/(pi*D*D/4)); % number of grains in each bin. this is simply geometry
% see equation 2 in the ripple paper from 1990


% Number of Impacts
dt = 1;
tmax = 1000;
t = 0:dt:tmax;
imax =length(t);


%Particle Impacts
for i=1:imax
    
particle_wind = 0.01*((xmax)*rand-x); %generate the particle trajectories

intercept = find(particle_wind < z); %find locations where the particle wind impacts the topo
%if (numel(intercept==0)) %in case sand is not hit

%else    
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

hold on
plot(x,particle_wind)
plot(x,z,'k')

%pause(0.1)
drawnow
%end
end