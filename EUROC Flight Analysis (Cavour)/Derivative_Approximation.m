%% Checking Convergence of the Method
clc; clear all;close all;

load('RocketPyData')

time = Rocketpydata(:,1);
velocity = Rocketpydata(:,3);
acceleration = Rocketpydata(:,4);

dVelocity = diff(velocity);  %calculating all the slopes of the curve, then plotting them together
dTime = diff(time);

%calculating derivative
acceleration_der = dVelocity./dTime;
acceleration_der_NOUDV =[];

for element = [1:size(acceleration_der)] %%putting acceleration = precedent value where NAN is calculated. It means same acceleration at same time
    if not(isnan(acceleration_der(element)))
        acceleration_der_NOUDV = [acceleration_der_NOUDV,acceleration_der(element)];
    else
        acceleration_der_NOUDV = [acceleration_der_NOUDV,acceleration_der(element+1)];
    end
end

acceleration_der_NOUDV= acceleration_der_NOUDV';


plot(time(2:end),acceleration_der_NOUDV,'r','LineWidth',1.5,'DisplayName',' Acceleration from Derivative')
hold on
plot(time,acceleration,'k--','LineWidth',1.5,'DisplayName',' Acceleration Simulated')
hold on
plot(time,velocity,'m','LineWidth',1.5,'DisplayName',' Acceleration Simulated')
hold off

%% Checking Convergence of the METHOD with EASYMINI DATA
clc; clear all;close all;

load("easyDATAwV.mat")
timeE = EASYdata(32:end,1);
accelerationE = EASYdata(32:end,2)
velocityE = EASYdata(32:end,4)

dv = diff(velocityE);
dt = diff(timeE);

acc = dv./dt;
acc_udv =[]
time_udv =[]


for element = [2:size(acc)] %%putting acceleration = 0 where NAN is calculated. It means same velocity at same time
    if not((isnan(acc(element))  || (acc(element)>200)))
        acc_udv = [acc_udv,acc(element)];
        time_udv =[time_udv,timeE(element)];
 
    end
end

plot(time_udv',acc_udv','r','LineWidth',1.5,'DisplayName',' Acceleration from Derivative')
hold on
plot(timeE,accelerationE,'k--','LineWidth',1.5,'DisplayName',' Acceleration by Easymini')
hold on
plot(timeE,velocityE,'b-.','LineWidth',1.5,'DisplayName',' Velocity Easymini')
hold off

legend("location","northeast")

%% Convergence of VEGA Data

clc
clear all
close all

load('easyDATAwV.mat')
load('RocketPyData')
load("vegaDATA")

time = Rocketpydata(:,1);
velocity = Rocketpydata(:,3);
acceleration = Rocketpydata(:,4);

timeE = EASYdata(32:end,1);
accelerationE = EASYdata(32:end,2)
velocityE = EASYdata(32:end,4)


%Acceleration From Vega
timeV = (VEGAdata(1:end,1)-VEGAdata(1,1))/10;  %removing offset, time in seconds
velocityV = VEGAdata(1:end,3);

%calculating f(x+1)-f(x) for each point
diffV = diff(velocityV); %this is a vector contaning all the differences [v(2:n,:) - v(1:n-1,:)].
diffT = diff(timeV); %vector containing all the differences [t(2:n,:) - t(1:n-1,:)]

acceleration_VEGA = diffV./diffT;
acceleration_VEGA_noUDV = []; 
time_VEGA_noUDV=[]

for element = [1:size(acceleration_VEGA)] %%putting acceleration = 0 where NON is showed. It means same velocity at same time
    if not(isnan(acceleration_VEGA(element)) || acceleration_VEGA(element)<-40)
        acceleration_VEGA_noUDV = [acceleration_VEGA_noUDV,acceleration_VEGA(element)];
        time_VEGA_noUDV = [time_VEGA_noUDV,timeV(element)]
        
    end
end
values =[]
for element = [1:size(acceleration_VEGA)]
    if acceleration_VEGA(element)==0
        values=[values,element];
    end
end

index = min(values);
display("apogee")
display(timeV(index))

% display("Maximum Acceleration VEGA");
% max(acceleration_VEGA_noUDV)

%Velocity
% plot(timeV,velocityV,'b','LineWidth',1.5,'DisplayName',' Velocity Vega')
% hold on
% plot(time,velocity,'g-.','LineWidth',1.5,'DisplayName',' Velocity simulation')
% hold on
% plot(timeE,velocityE,'m-.','LineWidth',1.5,'DisplayName',' Velocity easymini')

% %Acceleration
plot(time_VEGA_noUDV,acceleration_VEGA_noUDV,'r','LineWidth',1.5,'DisplayName',' Acceleration Vega')
hold on
plot(time,acceleration,'k--','LineWidth',1.5,'DisplayName',' Acceleration Simulated')
hold on
% plot(timeE,accelerationE,'b--','LineWidth',1.5,'DisplayName',' Acceleration Easymini')
% hold off

legend("location","northeast")
