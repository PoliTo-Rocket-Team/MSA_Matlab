% clear all
close all
clc

% WARNING
% this script needs mass_ts, thrust_ts and abcs_deploy to work
% generate them through step0_initialize_values.m or load rocketdata.mat

% INFO
% the scope of this file is simulating the trajectory of the rocket using
% the 'default' air brakes extension for the entirety of the air brakes
% control system activation time to calculate the apogee and the
% time-to-apogee used for trajectories generation in the
% step2_trajectories_sim_script.m script

% the apogee can be tuned by altering the mass of the rocket without the
% engine assembly (mass_struct) in the step0_initialize_values.m

% alternatively, the apogee_var.mat file can be loaded to load the required
% variables without the need of running this code


load MachExtension2Cdtable.mat
fprintf('MachExtension2Cd lookup table loaded \n\n');

default_abe = 0.5; % default air brakes extension

t_sim = 30; % simulation time (must be greater that time-to-apogee)

out = sim('simtime_calculator.slx');

apogee = max(out.vert_alt.Data);
t2a_index = find(out.vert_alt.Data == apogee);
t2a = out.vert_alt.Time(t2a_index);

mass_final = mass_ts.Data(end);

fprintf('Apogee is %.2f m at %.3f s (dry rocket mass %.3f kg) \n', apogee, t2a, mass_final);

target_apogee = 3000; % target apogee
margin = 0.001 * target_apogee; % 0.1% of the target apogee

if apogee >= (target_apogee - margin) && apogee <= (target_apogee + margin)
    fprintf('Apogee within 0.1%% margin of target (%.2f m), proceed \n\n', target_apogee);
else
    fprintf('Apogee outside 0.1%% margin of target (%.2f m), tune mass_struct in step0_initialize_values.m \n\n', target_apogee);
end

disp('t2a saved to workspace');
fprintf('apogee saved to workspace \n\n');


answer = input('do you want to save the variables to apogee_var.mat? (yes/[no]) \n->', 's');

if isempty(answer) || strcmpi(answer, 'no')
    fprintf('\n');
else
    save apogee_var.mat t2a apogee
    fprintf('\n');
    fprintf('t2a saved to apogee_var.mat \n');
    fprintf('apogee saved to apogee_var.mat \n\n');
end

clearvars -except thrust_ts mass_ts t y y_mat dy dy_mat surfaceBlockParameters1 t2a apogee abcs_deploy

disp('done, proceed to step 2');