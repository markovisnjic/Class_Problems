%% Surface Flow
%Marko Visnjic 8 March 2016
%% Parameters

clear all
figure (1)
clf

k = 10^-2;

g = 9.81; 

%% Topography
dx=2;
xmax=200;
x=0:dx:xmax;


zmax=20;
S = zmax/xmax;
z=zmax-(S*x);


%% Time
dt = 0.1; % in seconds
tmax = 3600; %in seconds
t = 0:dt:tmax;
nplots = 50;
tplot = tmax/nplots; % plotting interval

imax = length(t);

%% Precipitation
R0 = 0.005; 
period = 1800;
Ramp = 0.007;
precip = R0 + Ramp*sin(2*pi*t/period);
precip = precip/3600; % converts to m/second rates
precip = max(precip,0);

I0 = 1e-3; % infiltration if wet in m/hr
I0 = I0/3600; % converts to m/second
water = zeros(size(x));
h = water;
sub_zero = find(precip<0);
precip(sub_zero)=0;

zt = z+water; %total topo including water...spo water surface topo
ff = 0.5; % upslope weighting if>0.5...0.5 would be even weighting
vonkarman = 0.4;
z0 = 1e-4; % roughness constant in m

for i=1:imax
    wet = find(h>0);
    I=zeros(size(x));
    I(wet) = I0; % infiltration rate if wet
    P = precip(i)*ones(size(x)); %uniform rainfall
    hedge = ff*(h(1:end-1))+((1-ff)*h(2:end));
    dzdx = diff(z)/dx; %Slope
    
        ustar = sqrt(g*hedge.*abs(dzdx));
        ubar = zeros(size(hedge));
        flow = find(hedge>z0);
        ubar(flow) = (ustar(flow)./vonkarman).*(log(hedge(flow)/z0)-1); 
        q = ubar.*hedge;
    q = [0 q]; %No flux at top of hill
    dqdx = diff(q)/dx;
    dqdx = [dqdx dqdx(end)];

    dhdt = -dqdx + (P-I);
    h = h + (dhdt*dt);
    h = max(0,h); % eliminates negative water thicknesses
    zt = z + h;
    if(rem(t(i),tplot)==0)
        max(h);
    figure (1)
    subplot(3,1,1)
    plot(x,h*1000,'b')
    axis([0 xmax 0 5])
ylabel('Water Depth(mm)','FontSize',16,'fontname','arial')
        xlabel('Distance (m)','FontSize',16,'fontname','arial')
        drawnow
        hold off

    subplot(3,1,2)
    plot(x,q,'k')
        axis([0 xmax -2e-3 2e-3])
    ylabel('Q','FontSize',16,'fontname','arial')
        xlabel('Distance (m)','FontSize',16,'fontname','arial')
    hold off
    
    subplot(3,1,3)
    plot(t/3600,precip*3600*1000,'b') % precip in mm/hr
    axis([0 tmax/3600 0 1.1*max(precip)*3600*1000])
    ylabel('Time (hours)','FontSize',16,'fontname','arial')
    xlabel('Rainfall (mm/hr)','FontSize',16,'fontname','arial')

    end
    
end