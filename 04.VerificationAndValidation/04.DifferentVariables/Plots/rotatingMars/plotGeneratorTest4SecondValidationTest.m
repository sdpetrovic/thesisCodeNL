% This script lets me plot the 3-D graphs of the generated trajectories
% and plot the sensitivity analysis plots.
% For this the Cartesian generated data is used (also in the spherical
% computation case) and the latests files for the different variables.
%
% Stacha Petrovic 03-12-2016
% 
% version 1
% 
% Linux Ubuntu 16.04 LTS
%
%
% 
%
%

close all
clear all
clc

tic


%% Get the data

% Specify the paths
ThrustPath = fullfile('04.Version4SecondValidationTestThrust','ThrustAngles.csv'); % Create the path for the thrust angles file


ThrustAngles = csvread(ThrustPath); % Read the file



%% Generate the required vectors

% Provided thrust angles
Time = ThrustAngles(:,1); % Time [sec]
Altitude = ThrustAngles(:,2); % Altitude [km MOLA]
ThrustAzimuthAngles = ThrustAngles(:,3); % The thrust azimuth angles [deg]
ThrustElevationAngles = ThrustAngles(:,4); % The thrust elevation angles [deg]


% Assumed thrust angle model
fitTimeAzimuth = [0, 1, 3.09 6.81, 99.362, 99.363, 960];
fitThrustAzimuthAngles = [0, 0, 4.115, 4.154, 7.0835, 0, 0 ];
fitThrustAzimuthAnglesAltitude = [0, 0, 4.115, 4.154, 5.638, 7.0835, 0, 0 ];


fitTimeElevation = [0, 3.09, 9.965, 20.67, 97, 99.363, 960];
fitThrustElevationAngles = [0, 0, 0.174, 0.176, -21.45, 0, 0];
fitThrustElevationAnglesAltitude = [0, 0, 0.174, 0.176, -11.137, -21.45, 0, 0];


fitAltitudeAzimuth = [-0.6, -0.59, -0.5278, -0.277, 10.52 54.045, 54.046, 480];
fitAltitudeElevation = [-0.6, -0.5278, 0.0729, 2.176, 18.227, 51.5, 54.046, 480];





%% Plots


figure(1)
plot(Time,ThrustAzimuthAngles);
hold on
plot(Time,ThrustElevationAngles);
hold on
plot(fitTimeAzimuth,fitThrustAzimuthAngles);
hold on
plot(fitTimeElevation,fitThrustElevationAngles);

title('Thrust azimuth and elevation angles vs time'); % Give the figure a title
xlabel('Time [sec]'); % Label the different axes
ylabel('Thrust Azimuth/Elevation angle [deg]');

legend('Thrust Azimuth angle','Thrust Elevation angle','Fit Azimuth','Fit Elevation','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(2)
plot(Altitude,ThrustAzimuthAngles);
hold on
plot(Altitude,ThrustElevationAngles);
hold on
plot(fitAltitudeAzimuth,fitThrustAzimuthAnglesAltitude);
hold on
plot(fitAltitudeElevation,fitThrustElevationAnglesAltitude);

title('Thrust azimuth and elevation angles vs altitude'); % Give the figure a title
xlabel('Altitude [km MOLA]'); % Label the different axes
ylabel('Thrust Azimuth/Elevation angle [deg]');

legend('Thrust Azimuth angle','Thrust Elevation angle','Fit Azimuth','Fit Elevation','Location','NorthEastOutside'); % Add a legend in the top right corner



toc




