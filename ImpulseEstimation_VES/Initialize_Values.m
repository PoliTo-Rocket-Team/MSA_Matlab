clc
clearvars -except k* s t
close all
% thrust_table = readmatrix("AeroTech_O5500X-PS.csv");
% thrust = thrust_table(:,2);
% thrust_time = thrust_table(:,1);
thrust = s;
thrust_time = t;
time = linspace(0,28,1000);
thrust_interp = interp1(thrust_time,thrust,time); %interpolate where needed
thrust_interp(isnan(thrust_interp)) = 0;
% impulse_avg = 40960; %code to create a basic customized thrust curve
% thrust_avg = 8000;
% time = linspace(0,40,1000);
% t_burn = impulse_avg/thrust_avg;
% thrust = zeros(length(time),1);
% k = 1;
% while time(k)<t_burn
%     thrust(k) = thrust_avg;
%     k = k + 1;
% end
% delta = 0.75;
% rate = -thrust_avg/(2*delta);
% for k=1:length(time)
%     if (time(k)>t_burn-delta) && (time(k)<t_burn+delta)
%         thrust(k) = -rate*(t_burn+delta) + rate*time(k);
%     end
% end
thrust_ts = timeseries(thrust_interp,time);
mass_diff = 3.693;
mass_motor = 6.451;
t_burn = 3.6;
j = 1;
% for mass_struct = 15:2:30
        mass_struct = 19.23;
        m_init = mass_struct + mass_motor;
        m_final = m_init - mass_diff;
        k = -mass_diff/t_burn;
        i = 1;
        mass = ones(length(time),1)*m_final;
        while time(i) < t_burn
            mass(i) = k*time(i) + m_init;%linear decrease of total mass during burn time of motor
            i = i + 1;
        end
        mass_ts = timeseries(mass,time);
        out = sim("ThreeDOF_Rocket_Simulator_v2.slx");
        heights = getElement(out.yout,'Altitude');
        Velocity = getElement(out.yout,'Velocity');
        Acceleration = getElement(out.yout,'Acceleration');
        time_heights = heights.Values.Time;
        time_velocity = Velocity.Values.Time;
        time_acceleration = Acceleration.Values.Time;
        heights_val = heights.Values.Data - 1400; % 1400 initial altitude for SA Cup Launch
        Velocity_val = Velocity.Values.Data;
        Acceleration_val = Acceleration.Values.Data;
        km9_130_O5500(j) = max(heights_val);
        j = j+1;
%end

figure 
plot(time_heights,heights_val,"k","LineWidth",2)
xlabel("Time")
ylabel("Altitude")
grid on

figure 
plot(time_heights,Velocity_val,"k","LineWidth",2)
xlabel("Time")
ylabel("Velocity")
grid on

figure
plot(time_heights,Acceleration_val,"k","LineWidth",2)
xlabel("Time")
ylabel("Acceleration")
grid on
