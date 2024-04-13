%% startup
clc
 clear all
close all


%% MASS 
mass_table = readtable('Mass.txt');
time_mass = mass_table.Var1; %s
time_mass = [0;0.01;0.02;time_mass];
mass_vals = mass_table.Var2*1e-6; %kg
mass_vals = [12.86;12.86;12.86;mass_vals];
mass_ts = timeseries(mass_vals, time_mass);


%% INERTIA
inertia_table = readtable('inertia.txt');
time_inertia = inertia_table.Var1; %s
time_inertia = [0;time_inertia];
inertia_vals = inertia_table.Var2; %kg*m2
inertia_vals = [2.937;inertia_vals];
inertia_ts = timeseries(inertia_vals, time_inertia);

%% INERTIA TIME RATE
inertia_rate_vals = zeros(length(inertia_vals),1);
for i=[2:length(inertia_rate_vals)]
    h = time_inertia(i) - time_inertia(i-1);
    inertia_rate_vals(i-1) = (inertia_vals(i)-inertia_vals(i-1))/h; %kg*m2/s
end
inertia_rate_ts = timeseries(inertia_rate_vals, time_inertia);

%% THRUST
thrust_table = readtable('Thrust.txt');
time_thrust = thrust_table.Var1; %s
time_thrust = [0;time_thrust];
thrust_vals = thrust_table.Var2; %N
thrust_vals = [0;thrust_vals];
thrust_ts = timeseries(thrust_vals, time_thrust);

%% DRAG COEFFICIENT
cd_table = readtable('cd_time.txt');
time_cd = cd_table.Var1; %s
time_cd = [0;time_cd];
cd_vals = cd_table.Var2; 
cd_vals = [0.626;cd_vals];
cd_ts = timeseries(cd_vals, time_cd);

%% Static MARGIN
static = readtable('static_margin.txt');
time_st = static.Var1;
time_st = [0;time_st];
static_vals = static.Var2;
static_vals = [0;static_vals]
static_ts = timeseries(static_vals,time_st)



