clc
clearvars -except k*
close all
thrust_table = readmatrix("AeroTech_M2400T.csv"); %modified
thrust = thrust_table(:,2);
thrust_time = thrust_table(:,1);
time = linspace(0,27,1000); %modified
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
mass_motor = 6.451; %modified
motor_dry=2.758;
mass_diff=mass_motor-motor_dry;
t_burn = 3.3;
j = 1;
for mass_struct = 19.825 %
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
end
save("rocketdata.mat","thrust_ts","mass_ts")