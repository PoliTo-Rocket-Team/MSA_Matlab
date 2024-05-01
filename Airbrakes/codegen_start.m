clear all
close all
clc

% INFO
% this script automatically loads all the variables required for
% code generation without the need of generating them through the single
% scripts and opens the Simulink schematic of the code generator

load pid_var.mat
load rocketdata.mat
load trajectories.mat
load apogee_var.mat
load MachExtension2Cdtable.mat

fprintf('required workspace variables loaded \n\n');
disp('-> ready for code generation, opening codegen_model.slx... <-');

open codegen_model.slx
