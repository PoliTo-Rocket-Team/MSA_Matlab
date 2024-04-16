clear all
close all
clc

load rocketdata.mat
load trajectories.mat
load apogee_var.mat
load MachExtension2Cdtable.mat

disp('required workspace variables loaded');
disp('-> ready for simulation, opening airbrake_sim.slx... <-');

open airbrake_sim.slx
