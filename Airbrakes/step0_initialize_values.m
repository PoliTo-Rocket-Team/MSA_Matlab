clc
clear all
close all

% MASTER WARNING
% the PID is tuned for the parameters (if not previously modified by the
% user) which can be loaded through the sim_start.m script, and as such any
% modification in the variables which will be generated through this 3-step
% process might cause a suboptimal response

% therefore, in case of substantial modifications to the simulation
% variables, further tuning of the PID is required

% additional info on the PID tuning process can be found on the README.md
% located in the Airbrakes folder

% INFO
% the scope of this file is generating the 1x1 timeserie variables which
% correlate the thrust of the engine and the mass of the whole rocket
% assembly over time

% the thrust_ts and mass_ts timetables are needed to run
% step1_simtime_calculator_script.m, step2_trajectories_sim_script.m and
% ultimately the airbrake_sim.slx itself

% alternatively, the rocketdata.mat file can be loaded to load the
% aforementioned variables without the need of running this code

thrust_table = readmatrix("AeroTech_M2400T.csv"); % modified
thrust = thrust_table(:,2);
thrust_time = thrust_table(:,1);

time = linspace(0,30,1000); % modified
thrust_interp = interp1(thrust_time,thrust,time); % interpolate where needed
thrust_interp(isnan(thrust_interp)) = 0;

thrust_ts = timeseries(thrust_interp,time);
mass_motor = 6.451; % initial mass of the motor assembly
motor_dry=2.758; % final mass of the motor assembly
mass_diff=mass_motor-motor_dry;

t_burn = 3.3; % burn time of the engine
j = 1;

for mass_struct = 21.665 % mass of the rocket without the engine assembly
        m_init = mass_struct + mass_motor;
        m_final = m_init - mass_diff;
        k = -mass_diff/t_burn;
        i = 1;
        mass = ones(length(time),1)*m_final;
        while time(i) < t_burn
            mass(i) = k*time(i) + m_init; % linear decrease of total mass during burn time of motor
            i = i + 1;
        end
        mass_ts = timeseries(mass,time);
end


disp('thrust_ts saved to workspace');
fprintf('mass_ts saved to workspace \n\n');

answer = input('do you want to save the variables to rocketdata.mat? (yes/[no]) \n->', 's');

if isempty(answer) || strcmpi(answer, 'no')
    fprintf('\n');
else
    save rocketdata thrust_ts mass_ts
    fprintf('\n');
    fprintf('thrust_ts saved to rocketdata.mat \n');
    fprintf('mass_ts saved to rocketdata.mat \n\n');
end

clearvars -except thrust_ts mass_ts

disp('done, proceed to step 1');