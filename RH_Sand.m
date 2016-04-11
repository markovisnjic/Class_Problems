%%Grain bounce model
%% Rachel Havranek--- April,2016

%Initialize
D=0.25; %mm -- typical size of eolian sand
D=D/1000; %m
eta=0.35; %unitless
alpha=2;

xmax=10; %10 m array size
%dx=5*D;
dx=.02;
x=0:dx:xmax;

tmax=10; %days
tmax=tmax*24*60*60; %seconds
dt=tmax/100000;
t=0:dt:tmax;
imax=length(t);
nplots = 100; %so our code runs faster - we'll plot one hundred times
tplot = tmax/nplots; %so our code runs faster - divides total time up into 100 time steps

a=19999;
b=20001;
Nint=round(a + (b-a).*rand(1,(length(x))));
N=Nint.*(ones(size(x)));

%zintercept=.5+rand;
nframe = 0;

%%RUN
for i=1:imax
   %height of each bin
   z=(pi*N*D^2)/(4*(1-eta)*dx);
   
   %impacting grains
%    c=7;
%    d=10;
%    alpha=c+(d-c).*rand(1,1);
min_intercept=z(1);
z_intercept=x*tand(alpha)+z;
max_intercept=max(z_intercept);
intercept_range=max_intercept-min_intercept;

intercept=rand*intercept_range +min_intercept;
   ximpact=intercept-x*(tand(alpha));
   
   %find the index for this value%
   potential_impacts=find(z>ximpact);
   impact=potential_impacts(1);
   
   %number of grains moved
   grm=(1+(round(20*rand)));
    
  N(impact)=N(impact)-grm;

   
  if impact>=(501-25)
      impact=(impact-500);
  end
  N(impact+25)=N(impact+25)+grm;
   
   if (rem(t(i),tplot)==0)
    nframe = nframe+1;
   figure(1)
plot(x,z,'k','linewidth',1.5)
axis ([0,xmax,0,.1])
    pause(0.5)
   

   end
end
