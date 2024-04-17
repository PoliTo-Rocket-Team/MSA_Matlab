clear all
close all
clc

% INFO
% this script automatically loads all the variables required for
% simulation without the need of generating them through the single
% scripts, opens the Simulink schematic of the air brakes control system
% and runs it once

load pid_var.mat
load rocketdata.mat
load trajectories.mat
load apogee_var.mat
load MachExtension2Cdtable.mat

fprintf('required workspace variables loaded \n\n');
disp('-> ready for simulation, running airbrake_sim.slx... <-');

open Airbrake_sim.slx
out = sim('Airbrake_sim.slx');