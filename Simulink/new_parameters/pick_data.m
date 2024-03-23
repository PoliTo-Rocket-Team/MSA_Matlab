
clc;
clear all;
close all;
load("flightData.csv")
flightData = flightData()
%% Mass
fileID = fopen('mass.txt', 'w');
if fileID == -1
    error('Unable to open file for writing');
end
mass = [flightData(1:510,1),flightData(1:510,2).*1000]
% Write the data to the file
fprintf(fileID, '%f %f\n', mass');
% Close the file
fclose(fileID);
%% Thrust
fileID = fopen('thrust.txt', 'w');
if fileID == -1
    error('Unable to open file for writing');
end
thrust = [flightData(1:510,1),flightData(1:510,7)]
% Write the data to the file
fprintf(fileID, '%f %f\n', thrust');
% Close the file
fclose(fileID);

%%  Drag Coefficient
%with time
fileID = fopen('cd_time.txt', 'w');
if fileID == -1
    error('Unable to open file for writing');
end
thrust = [flightData(1:510,1),flightData(1:510,8)]
% Write the data to the file
fprintf(fileID, '%f %f\n', thrust');
% Close the file
fclose(fileID);

%% Drag coefficient
%%with mach
fileID = fopen('cd_mach.txt', 'w');
if fileID == -1
    error('Unable to open file for writing');
end
drag = [flightData(1:510,6),flightData(1:510,8)]
dlmwrite('cd_mach.txt', drag, 'delimiter', ',');
% Close the file
fclose(fileID);

%% Inertia
fileID = fopen('inertia.txt', 'w');
if fileID == -1
    error('Unable to open file for writing');
end
inertia = [flightData(1:510,1),flightData(1:510,3)]
% Write the data to the file
fprintf(fileID, '%f %f\n', inertia');
% Close the file
fclose(fileID);

%% STATIC MARGIN
fileID = fopen('static_margin.txt', 'w');
if fileID == -1
    error('Unable to open file for writing');
end
inertia = [flightData(1:510,1),flightData(1:510,3)*0.000254] %seconds and meters
% Write the data to the file
fprintf(fileID, '%f %f\n', inertia');
% Close the file
fclose(fileID);