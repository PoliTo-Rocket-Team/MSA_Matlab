clc
clearvars -except O*
close all
thrust_table = readmatrix("Cesaroni_40960O8000-P.csv");
thrust = thrust_table(:,2);
thrust_time = thrust_table(:,1);
time = linspace(0,40,1000);
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
mass_diff = 18.6;
t_burn = 5.1;
% mass_struct = 21.9;
O8000 = zeros(8,1);
j = 1;
for m_init = 45:5:80
    if m_init-mass_diff > 0
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
        heights_val = heights.Values.Data;
        O8000(j) = max(heights_val);
    else
        O8000(j) = 0;
    end
    j = j + 1;
end