% clear all
% close all
% clc

%using curvefit to export as a simulink lookup table
%table parameter is saved MachExtension2Cdtable.mat

%% Extension as a function of Cd and Mach
T=readtable("Rocketteam_CD_MD.xlsx"); %from missile datcom
T_mat=table2array(T);
Extension=T_mat(:,1);
CD=T_mat(:,2);
Mach=T_mat(:,3);

[xData, yData, zData] = prepareSurfaceData( Mach, CD, Extension );

ft = 'cubicinterp'; %cubic spline interpolation
opts = fitoptions( 'Method', 'CubicSplineInterpolant' );
opts.ExtrapolationMethod = 'none';
opts.Normalize = 'on';

% extensionfit = fit([xData,yData],zData,'poly33','normalize','on');
[extensionfit, gof] = fit( [xData, yData], zData, ft, opts );
figure(1)
plot(extensionfit,[xData,yData],zData);% xlabel( 'Mach', 'Interpreter', 'none' );
xlabel( 'Mach', 'Interpreter', 'none' );
ylabel( 'Cd', 'Interpreter', 'none' );
zlabel( 'Extension', 'Interpreter', 'none' );
grid on















