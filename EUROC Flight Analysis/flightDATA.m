
%% VELOCITY AND ALTITUDE VEGA
clc; 
clear all;
close all;

load("vegaDATA")
load("easyDATAwV")

%VEGA
timeV = (VEGAdata(1:end,1))/10  %removing offset, time in seconds
altitudeV = VEGAdata(1:end,2)
velocityV = VEGAdata(1:end,3)

%EASYMINI
timeE = EASYdata(32:end,1);
accelerationE = EASYdata(32:end,2)
altitudeE = EASYdata(32:end,3)
velocityE = EASYdata(32:end,4)
voltageD = EASYdata(32:end,5)
voltageM = EASYdata(32:end,6)

figure(1)
yyaxis left;
plot(timeV,altitudeV,'b','DisplayName','altitude VEGA')
hold on
ylabel("Altitude(m)")
ylim([-70,3000])

yyaxis right
plot(timeV,velocityV,'r','DisplayName','velocity VEGA')
hold off
ylabel("Velocity (m/s)")
ylim([-90,350])

xlabel('time(S)') 
grid on
legend('Location','northeast')

%% ALTITUDE COMPARISON VEGA EASYMINI and SIMULATION
clc; clear all; close all;

clc; 
clear all;
close all;

load("vegaDATA")
load("easyDATAwV")
load("RocketPyData")

%VEGA
timeV = (VEGAdata(2:end,1)-VEGAdata(2,1))/10;  %removing offset, time in seconds
altitudeV = VEGAdata(2:end,2);
display("apogee vega")
max(altitudeV)

%EASYMINI
timeE = EASYdata(32:end,1);
altitudeE = EASYdata(32:end,3);
display("apogee easymini")
max(altitudeE)

%ROCKETPY
time_s = Rocketpydata(2:end,1);
ground_level =  Rocketpydata(2,2)
altitude_s = Rocketpydata(2:end,2)-ground_level;

display("apogee RocketPy")
max(altitude_s)

ylabel("Altitude(m)")
xlabel("Time(s)")

plot(timeE,altitudeE,'b','DisplayName','altitude Easymini');
hold on
plot(timeV,altitudeV,'r','DisplayName','altitude Vega');
hold on
plot(time_s,altitude_s,'k','DisplayName','altitude RocketPy');
grid on
legend('Location','northeast')

%% VELOCITY COMPARISON BETWEEN VEGA EASYMINI AND ROCKETPY
clc; 
clear all;
close all;

load("vegaDATA")
load("easyDATAwV")
load("RocketPyData")

%VEGA
timeV = (VEGAdata(1:end,1)-VEGAdata(1,1))/10;  %removing offset, time in seconds
velocityV = VEGAdata(1:end,3);



%EASYMINI
timeE = EASYdata(32:end,1);
velocityE = EASYdata(32:end,4);

%ROCKETPY
time_s = Rocketpydata(2:end,1);
velocity_s = Rocketpydata(2:end,3);

display("velocity simulated with RocketPy")
max(velocity_s)

display("max velocity easymini")
max(velocityE)



plot(timeE,velocityE,'b','LineWidth',1.5,'DisplayName','velocity Easymini')
hold on
plot(timeV,velocityV,'k','LineWidth',1.5,'DisplayName','velocity Vega')
hold on
plot(time_s,velocity_s,'r','LineWidth',1.5,'DisplayName','velocity simulated by RocketPy')
hold off
grid on
legend('Location','northeast')

ylabel("velocity (m/s)")
xlabel("Time(s)")




%% ACCELERATION COMPARISON BETWEEN EASYMINI AND ROCKETPY
clc; 
clear all;
close all;

load("easyDATAwV")
load("RocketPyData")

%EASYMINI
timeE = EASYdata(32:end,1);
accelerationE = EASYdata(32:end,2);

%ROCKETPY
time_s = Rocketpydata(2:end,1);
acceleration_s = Rocketpydata(2:end,4);

display("max acceleration rocketpy");
max(acceleration_s/9.81)
display("max acceleration easymini");
max(accelerationE)

relative_error = (abs(max(acceleration_s)-max(accelerationE))/max(acceleration_s))*100
plot(timeE,accelerationE,'b','LineWidth',1.5,'DisplayName','acceleration Easymini')
hold on
plot(time_s,acceleration_s,'r--','LineWidth',1.5,'DisplayName','acceleration simulated by RocketPy')
hold off

grid on
legend('Location','northeast')

ylabel("Acceleration (m/s^2)")
xlabel("Time(s)")

%% MAIN AND APOGEE VOLTAGES EASYMINI
clc; clear all; close all;

load("vegaDATA")
load("easyDATAwV")

%VEGA
timeV = (VEGAdata(2:end,1)-VEGAdata(2,1))/10  %removing offset, time in seconds
altitudeV = VEGAdata(2:end,2)-VEGAdata(2,2)
velocityV = VEGAdata(2:end,3)-VEGAdata(2,3)

%EASYMINI
timeE = EASYdata(32:end,1);
accelerationE = EASYdata(32:end,2)
altitudeE = EASYdata(32:end,3)-EASYdata(32,3)
velocityE = EASYdata(32:end,4)
voltageD = EASYdata(32:end,5)
voltageM = EASYdata(32:end,6)


plot(timeE,voltageD,'b','DisplayName','Drogue Voltage')
hold on
plot(timeE,voltageM,'r','DisplayName','Main Voltage')
hold off
grid on
legend('Location','northeast')

ylabel("voltage(V)")
xlabel("Time(s)")


%% MAIN AND APOGEE FLIGHT STATUS CATS VEGA
clc; clear all; close all;

clc; 
clear all;
close all;

load("vegaDATA")
load("easyDATAwV")

%VEGA
timeV = (VEGAdata(1:end,1)-VEGAdata(1,1))/10  %removing offset, time in seconds
altitudeV = VEGAdata(1:end,2)
velocityV = VEGAdata(1:end,3)
drogue = VEGAdata(1:end,4)
main = VEGAdata(1:end,5)

%EASYMINI
timeE = EASYdata(32:end,1);
accelerationE = EASYdata(32:end,2)
altitudeE = EASYdata(32:end,3)-EASYdata(32,3)
velocityE = EASYdata(32:end,4)
voltageD = EASYdata(32:end,5)
voltageM = EASYdata(32:end,6)

plot(timeV,drogue,'b','LineWidth',1.5,'DisplayName','Drogue Status')
hold on
plot(timeV,main,'r','LineWidth',1.5,'DisplayName','Main Status')
hold on
plot(timeE,voltageD,'g','DisplayName','Drogue Voltage')
hold on
plot(timeE,voltageM,'k','DisplayName','Main Voltage')
hold off


grid on
legend('Location','southwest')

ylabel("status(0/1) or Voltage(V)")
xlabel("Time(s)")

%% Calculation of DRAG FORCE 
clc; clear all; close all;

load("easyDATAwV")
load("Boost_Data")
load("acceleration_SIM")

%Mass of the rocket

%dry_mass
m_0 = 8.219
%wet mass
m = 10.132

%COASTING
timeE = EASYdata(404:2309,1);
accelerationE = EASYdata(404:2309,2)
time_SIM = Acceleration_SIM(1:324,1)
acceleration_SIM = Acceleration_SIM(1:324,2)*9.81

%BURNOUT
time = Boost_DATA(:,1)
thrust = Boost_DATA(:,2)
acceleration = Boost_DATA(:,3)
acceleration_sim = Boost_DATA(:,4).*9.81

%BURNOUT
f_mass = (m-m_0)/3.8
mt = @(t) m - f_mass*t  %where t can be maximum 3.8s, time of burnout

%Real Boost Phase
Drag_during_burnout =@(t)  thrust- mt(t).*(acceleration+9.81)
%Simulated Boost Phase
Drag_during_burnout_SIM = @(t)  thrust- mt(t).*(acceleration_sim+9.81)

%COASTING

%Real Drag Coasting Phase
Drag_coasting = - m_0*(accelerationE+9.81)
%Simulated Drag Coasting Phase
Drag_coasting_SIM = - m_0*(acceleration_SIM+9.81)

%INITIALIZING DRAG VECTORS
DragVector = Drag_during_burnout(time) 
DragVector_SIM = Drag_during_burnout_SIM(time) 


DragVector =[DragVector',Drag_coasting']'   %adding to DragVector the drag coasting
DragVector_SIM =[DragVector_SIM',Drag_coasting_SIM']'



%Setting Drag to 0 when  is less than 0
for element = [1:size(DragVector)]
    if DragVector(element)<0
        DragVector(element)=0;
    end
end

for element = [1:size(DragVector_SIM)]
    if DragVector_SIM(element)<0
        DragVector_SIM(element)=0;
    end
end


timetot = [time',timeE']'  %combining time
timetot_SIM = [time',time_SIM']'

figure(1)

plot(timetot,DragVector,'k','LineWidth',1.5,'DisplayName',' Real Drag Force during boostphase and coasting')
hold on
plot(timetot_SIM,DragVector_SIM','b','LineWidth',1.5,'DisplayName','Simulated Drag Force ')

grid minor
legend("location","northeast")



%% Comparison Acceleration Simulated and Real Acceleration Easymini

clc; clear all; close all;

load("easyDATAwV")
load("Boost_Data")
load("acceleration_SIM")

%COASTING
timeE = EASYdata(404:2309,1);
accelerationE = EASYdata(404:2309,2)
time_SIM = Acceleration_SIM(1:324,1)
acceleration_SIM = Acceleration_SIM(1:324,2)*9.81

%BURNOUT
time = Boost_DATA(:,1);
acceleration = Boost_DATA(:,3);
acceleration_sim = Boost_DATA(:,4).*9.81;

acceleration = [acceleration',accelerationE']'
acceleration_sim = [acceleration_sim',acceleration_SIM']'
display("Maximum Acceleration Simulated")
max(acceleration_sim)
time_sim = [time',time_SIM']'
time_Real = [time',timeE']'

figure(1)
plot(time_sim,acceleration_sim,'k','LineWidth',1.5,'DisplayName','Simulated acceleration from Burnout to Coasting')
hold on
plot(time_Real',acceleration','b','LineWidth',1.5,'DisplayName','Real acceleration from Burnout to Coasting')
legend("location","northeast")
hold off
grid on


%% DRAG FORCE CALCULATION
clc
clear all
close all

load("OpenRocket_BoostData.mat")
load('easyDATAwV.mat')
load('RocketPyData')

%dry_mass
m_0 = 8.219;
%wet mass
m = 10.132;

%BURNOUT
time_EM = (linspace(0,3.80,380))';
acc_EM = EASYdata(2:381,2);
thrust_s = spline(time,thrust,time_EM);
thrust_s = [thrust_s(8:end) ; zeros(7,1)];

%COASTING
timeE = EASYdata(404:2309,1);
accelerationE = EASYdata(404:2309,2)

f_mass = (m-m_0)/3.8;
mt = @(t) m - f_mass*t ;
Drag_during_burnout =@(t)  thrust_s - mt(t).*(acc_EM+9.81);
DragVector = Drag_during_burnout(time_EM);  %during burnout

Drag_coasting = - m_0*(accelerationE+9.81)

time_s = Rocketpydata(2:349,1)
drag_s = Rocketpydata(2:349,5)

DragVector =[DragVector',Drag_coasting']'   %adding to DragVector the drag coasting
timeE =[timeE',time_EM']'

for element = [1:size(DragVector)]
    if DragVector(element)<0
        DragVector(element)=0;
    end
end

for element = [1:size(drag_s)]
    if drag_s(element)<0
        drag_s(element)=0;
    end
end

plot(timeE,DragVector,'k','LineWidth',1.5,'DisplayName',' Real Drag Force during boostphase and coasting')
hold on 
plot(time_s,drag_s,'r','LineWidth',1.5,'DisplayName','  Drag Force simulated by RocketPy')
hold off

legend("location","northeast")
