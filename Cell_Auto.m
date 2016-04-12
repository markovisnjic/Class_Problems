%%Cellular Automation random walk
%Particle walks around a box
%Marko Visnjic
figure(1)
clf
%Generates a Random walking cell
xmax = 50;
x = 0:xmax;
ymax = 50;
y = 0:ymax;

%Initial Surface
surface = zeros([xmax  ymax]);
particle = 1; % Particle is a "1" empty space is "0"
Init_PosX= xmax/2; %Starting Position
Init_PosY= ymax/2; %Starting Position
surface(Init_PosX,Init_PosY)= particle; %Starting the particle

%%Time
tmax = 100;
t = 0:tmax;
imax = length(t);



for i = 1:imax
    move = rand; %Every iteration a random number will be generated to correspond with moving the particle around
    if move <= 0.25; %Move up one in X
        motionX = 1;
        motionY = 0;
        
    elseif (0.25 < move) && (move <= 0.5) %Move down one in X
            motionX = -1;
            motionY = 0;
    end
    if (0.5 < move) && (move <= 0.75) %Move down one in Y
            motionY = -1;
            motionX = 0;
            
    elseif  move > 0.75 %Move down one in X
            motionY = 1;
            motionX = 0;
    end
    
    %Update the particle position
    Init_PosX = Init_PosX + motionX; 
    Init_PosY = Init_PosY + motionY; 
    
    if Init_PosX >= xmax
     Init_PosX = xmax;
    end
    if Init_PosY >= ymax
     Init_PosY = ymax;
    end
    
    %Set the new position equal to the particle  
    surface(Init_PosX+motionX,Init_PosY+motionY)= particle;
  
    image(surface,'CDataMapping','scaled')
    %Reset the surface and start again
  surface(:,:)=0;
  pause(0.1)
end
