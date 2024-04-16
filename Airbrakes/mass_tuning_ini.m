close all
clear all
clc

%% Cd function
T_Cd=readtable("Rocketteam_CD_MD.xlsx");
T_Cd_mat=table2array(T_Cd);
i=1; %select 50% extension
Cd=flip(T_Cd_mat(1+20*(i-1):20+20*(i-1),2));
Mach=flip(T_Cd_mat(1+20*(i-1):20+20*(i-1),3));

%%
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
m_init=zeros(100,1);
mass_struct = linspace(19,21,100);
for j=1:length(mass_struct) %
        m_init(j) = mass_struct(j) + mass_motor;
        m_final(j) = m_init(j) - mass_diff;
        k = -mass_diff/t_burn;
        i = 1;
        mass{j} = ones(length(time),1)*m_final(j);
        while time(i) < t_burn
            mass{j}(i) = k*time(i) + m_init(j);%linear decrease of total mass during burn time of motor
            i = i + 1;
        end
        mass_ts{j} = timeseries(mass{j},time);
        out = sim("mass_tuning.slx");
        heights = getElement(out.yout,'Altitude');
        heights_val{j} = heights.Values.Data;
        km3_M2400T(j) = max(heights_val{j}); %modified
end
mass_struct(53)