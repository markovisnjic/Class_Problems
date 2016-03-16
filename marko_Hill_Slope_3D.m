%% Hill Slope 3D
% Written by Marko Visnjic 15/3/16

figure (2)
clf

figure(3)
clf

%% Parameters
kappa = 0.0175; %Topographic Diffusivity m^2/year
k = 29.75;
rho_r = 2700; %Density of Rock
rho_b = 1700; %Bulk Density of Regolith
A0 = 50*10^-6; %Weathering of Bare Rock
A1 = 50*10^-6; %Regolith Dependent
b = 2*10^-4; %Thickness; no, b is not a thickness...look at its dimensions
h1 = 0.5; %Charectaristic Height; 
%characteristic soil depth

%% Arrays 

% X component

dx = 1;
xmax = 100;
x = -xmax:dx:xmax;
xmid = x(1:end-1)+(dx/2);

xone = ones(size(x));
% Y component

dy = 1;
ymax = 80; % making it look different so that we can see which axis is which
y = -ymax:dy:ymax;
ymid = y(1:end-1)+(dy/2);

yone = ones(size(y));
% Z component

dz = 0.5;
zmax = 200;
z = 0:dz:zmax;
slope = 0.1;

%zmid = z(1:end-1)+(dz/2);
%Bedrock Topo

[X,Y] = meshgrid(x,y);
[Xmid,Ymid] = meshgrid(xmid,ymid);

Sx = 0.01;
Sy = 0.01;
zb = zmax - sqrt(Sx*(X.^2)+Sy*(Y.^2));
zb0 = zb; %records initial conditions

% test to make sure your weathering soil production algorithm works
htest = 0:0.01:1;
wtest = min(A0+(b.*htest),A1*exp(-htest/h1)); % Weathering Rate
figure(3)
plot(htest,wtest)


%Soil Array
%[Xone,Yone] = meshgrid(xone,yone);
h0 = 0.01;
%h = h0.*Xone;
h = h0*ones(size(X));

Z = zb+h;
%Time
nplot = 50;
dt = 20;
tmax = 10000;
t = 0:dt:tmax;
tplot = tmax/nplot;
imax = length(t);

ci = (-0.25e-4)*ones(size(t));
%% Solve Diffusion Eq.
for i = 1:imax
    dzdx = diff(Z,1,2)/dx;
    dzdy = diff(Z,1,1)/dy;
   
    w = min(A0+(b.*h),A1*exp(-h/h1)); % Weathering Rate
    while size(w)<[162 202]
    wpad_1 = w(1,:);
    w = vertcat(wpad_1,w);
    wpad_2 = w(:,1);
    w = horzcat(wpad_2,w);
    end
    
    %Flux X
    qx = -k*dzdx; 
    
    qxpad_1 = qx(1,:); 
    qxpad_2 = vertcat(qx(1),qx(:,end));
    qx = vertcat(qxpad_1, qx);
    qx = horzcat(qx,qxpad_2);
    
    %Flux Y
    qy = -k*dzdy; 
    %Padding Flux
    qypad_1 = qy(1,:); 
    qypad_2 = vertcat(qy(1),qy(:,end));
    qy = vertcat(qypad_1, qy); 
    qy = horzcat(qy, qypad_2);
    
    dqxdx = diff(qx)/dx;
    
    dqydy = diff(qy)/dy;
    
    dQxdx = diff(qx,1,2)/dx;
    Qxpad_1 = dQxdx(:,end);
    Qxpad_2 = dQxdx(:,1);
    dQxdx = horzcat(Qxpad_2,dQxdx,Qxpad_1);
    
    dQydy = diff(qy,1,1)/dy;
    Qypad_1= dQxdx(end,:);
    Qypad_2 = dQxdx(1,:);
    dQydy = vertcat(Qypad_2, dQydy,Qypad_1);
    
    dhdt = w*(rho_r/rho_b)-(1/rho_b)*dQxdx-(1/rho_b)*dQydy;
  while size(h)<size(w)
    hbc_2 = h(:,end);
    h = [h hbc_2];
    hbc_1 = h(1,:);
    h = vertcat(hbc_1,h);
    h = h + (dhdt*dt);
  end  
    h = max(h,0);
    
   %Boundary Conditions
   
    while size(zb)<size(w)
    zbpad_1 = zb(:,1);
    zb = [zbpad_1 zb];
    zbpad_2 = zb(end,:);
    zb = vertcat(zbpad_2,zb);
    end
    
    zb(:,2:end-1) = zb(:,2:end-1) - (w(:,2:end-1)*dt);
    zb(:,1) =  zb(:,1) - ci(i)*dt;
    zb(:,end) = zb(:,end) - ci(i)*dt;
    
    zb(2:end-1,:) = zb(2:end-1,:) - (w(2:end-1,:)*dt);
    zb(1,:) =  zb(1,:) - ci(i)*dt;
    zb(end,:) = zb(end,:) - ci(i)*dt;
    
    
    h(:,1) = 0;
    h(:,end) = 0;
    h(end,:)=0;
    h(1,:)=0;
    
    z = zb+h;
    
    while size(X)<size(w)
    Xpad_1 = X(1,:);
    X=vertcat(Xpad_1,X);
    Xpad_2 = X(:,end);
    X = horzcat(Xpad_2,X);
    end
    while size(Y)<size(w)
    Ypad_1 = Y(1,:);
    Y=vertcat(Ypad_1,Y);
    Ypad_2 = Y(:,end);
    Y = horzcat(Ypad_2,Y);
    end
  if (rem(t(i),tplot)==0)
    figure(2)
    surf(X,Y,z)
    hold on
    surf(X,Y,zb)
    set(gca,'XDIR','reverse')
    view(-210,30)
    axis([-xmax xmax -ymax ymax min(min(zb0)) zmax])
    colormap(parula)
    colorbar('eastoutside')
    shading interp
    hidden off
    set(gca,'fontsize',18,'fontname','arial');
    xlabel('Distance_X[m]','fontsize',18,'fontname','arial');
    ylabel('Distance_Y[m]','fontname','arial','fontsize',21) ;
	zlabel('Height','fontname','arial','fontsize',21);
    drawnow
    hold off
    end
end


