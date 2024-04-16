clear all
close all
clc

load rocketdata.mat
load trajectories.mat
load apogee_var.mat
load MachExtension2Cdtable.mat

fprintf('required workspace variables loaded \n\n');
disp('-> ready for simulation, opening airbrake_sim.slx... <-');

open airbrake_sim.slx
