%%Sand Ripple Evolution 1D
%Marko Visnjic

clear all
figure(1)

%% Arrays

D = 0.005; % grain diameter
dx = 0.05; % make these bins much bigger than your grains
xmax = 2;
x = 0:dx:xmax;

% here is how i propose you make your initial topo and 
% back out how many grains are in each bin
z0=0.1;
z = z0+0.01*rand(size(x)); % random topography with a st dev 0.01

porosity = 0.35;
N = ceil(dx*z*(1-porosity)/(pi*D*D/4)); % number of grains in each bin. this is simply geometry
% see equation 2 in the ripple paper from 1990


% Time (actually, numer of impacts)
dt = 1;
tmax = 1000;
t = 0:dt:tmax;
imax =length(t);

%Initial Topography
%zb = x;

%X bins

%Particle Impacts
for i=1:imax
particle_wind = rand(z)-x; %generate the particle trajectories