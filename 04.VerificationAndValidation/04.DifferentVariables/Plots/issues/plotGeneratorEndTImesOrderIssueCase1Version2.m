% This script lets me plot the 3-D graphs of the generated trajectories
% and plot the sensitivity analysis plots.
% For this the Cartesian generated data is used (also in the spherical
% computation case) and the latests files for the different variables.
%
% Stacha Petrovic 18-12-2016
% 
% version 2
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

format long

tic


%% Get the data

% pathToValidationFolder = '/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/04.VerificationAndValidation/04.DifferentVariables';
% 
% currentVariableFolder = 'Order';
% currentVariableUnit = '[-]';
% currentCaseFolder = 'Case1/issue';


% Automatic latest file name (Source:
% https://www.youtube.com/watch?v=D8UkAOhsyDI)

% Write the paths to the folders
% pathTSI = fullfile(pathToValidationFolder,'TSI',currentVariableFolder,currentCaseFolder);
pathTSI = fullfile(pwd,'/EndTimes');
% pathRKF = fullfile(pathToValidationFolder,'RKF',currentVariableFolder,currentCaseFolder);

% Find all files ending with .csv in that folder
dTSI = dir(fullfile(pathTSI,'*.csv'));
% dRKF = dir(fullfile(pathRKF,'*.csv'));

% Find the latest file
datesTSI = [dTSI.datenum]; % Get all the dates
[~,newestIndexTSI] = max(datesTSI); % Get the index of the newest date
newestFileTSIOrder20 = dTSI(newestIndexTSI-2); % Get the data corresponding to the newest file
newestFileNameTSIOrder20 = newestFileTSIOrder20.name; % Get the name of the newest file

newestFileTSIOrder30 = dTSI(newestIndexTSI); % Get the data corresponding to the newest file
newestFileNameTSIOrder30 = newestFileTSIOrder30.name; % Get the name of the newest file

newestFileTSIOrder29 = dTSI(newestIndexTSI-1); % Get the data corresponding to the newest file
newestFileNameTSIOrder29 = newestFileTSIOrder29.name; % Get the name of the newest file

% 
% datesRKF = [dRKF.datenum]; % Get all the dates
% [~,newestIndexRKF] = max(datesRKF); % Get the index of the newest date
% newestFileRKF = dRKF(newestIndexRKF); % Get the data corresponding to the newest file
% newestFileNameRKF = newestFileRKF.name; % Get the name of the newest file


% dTSI.date




% Specify the paths
VariablePathTSIOrder20 = fullfile(pathTSI,newestFileNameTSIOrder20); % Create the path for the TSI data
VariablePathTSIOrder30 = fullfile(pathTSI,newestFileNameTSIOrder30); % Create the path for the TSI data
VariablePathTSIOrder29 = fullfile(pathTSI,newestFileNameTSIOrder29); % Create the path for the TSI data

% VariablePathRKF = fullfile(pathRKF,newestFileNameRKF); % Create the path for the RKF data


VariableVectorTSIOrder20 = csvread(VariablePathTSIOrder20); % Read the file
VariableVectorTSIOrder30 = csvread(VariablePathTSIOrder30); % Read the file
VariableVectorTSIOrder29 = csvread(VariablePathTSIOrder29); % Read the file

% VariableVectorRKF = csvread(VariablePathRKF); % Read the file







%% Generate the required vectors

% Desired data
timeTSIOrder20 = VariableVectorTSIOrder20(:,6); % Time


% States
xPositionTSIOrder20 = VariableVectorTSIOrder20(:,7); % x-position TSI [km]
yPositionTSIOrder20 = VariableVectorTSIOrder20(:,8); % y-position TSI [km]
zPositionTSIOrder20 = VariableVectorTSIOrder20(:,9); % z-position TSI [km]
xVelocityTSIOrder20 = VariableVectorTSIOrder20(:,10); % x-velocity TSI [km/s]
yVelocityTSIOrder20 = VariableVectorTSIOrder20(:,11); % y-velocity TSI [km/s]
zVelocityTSIOrder20 = VariableVectorTSIOrder20(:,12); % z-velocity TSI [km/s]
MassTSIOrder20 = VariableVectorTSIOrder20(:,13); % Mass TSI [kg]

% Desired data
timeTSIOrder30 = VariableVectorTSIOrder30(:,6); % Time


% States
xPositionTSIOrder30 = VariableVectorTSIOrder30(:,7); % x-position TSI [km]
yPositionTSIOrder30 = VariableVectorTSIOrder30(:,8); % y-position TSI [km]
zPositionTSIOrder30 = VariableVectorTSIOrder30(:,9); % z-position TSI [km]
xVelocityTSIOrder30 = VariableVectorTSIOrder30(:,10); % x-velocity TSI [km/s]
yVelocityTSIOrder30 = VariableVectorTSIOrder30(:,11); % y-velocity TSI [km/s]
zVelocityTSIOrder30 = VariableVectorTSIOrder30(:,12); % z-velocity TSI [km/s]
MassTSIOrder30 = VariableVectorTSIOrder30(:,13); % Mass TSI [kg]

% Desired data
timeTSIOrder29 = VariableVectorTSIOrder29(:,6); % Time


% States
xPositionTSIOrder29 = VariableVectorTSIOrder29(:,7); % x-position TSI [km]
yPositionTSIOrder29 = VariableVectorTSIOrder29(:,8); % y-position TSI [km]
zPositionTSIOrder29 = VariableVectorTSIOrder29(:,9); % z-position TSI [km]
xVelocityTSIOrder29 = VariableVectorTSIOrder29(:,10); % x-velocity TSI [km/s]
yVelocityTSIOrder29 = VariableVectorTSIOrder29(:,11); % y-velocity TSI [km/s]
zVelocityTSIOrder29 = VariableVectorTSIOrder29(:,12); % z-velocity TSI [km/s]
MassTSIOrder29 = VariableVectorTSIOrder29(:,13); % Mass TSI [kg]


% 
% 
% % Difference in metres w.r.t. the order 20 end state
% 
% xPositionTSIdifference = abs(TruthTSI(1)-xPositionTSI)*1000; % x-position TSI [m]
% yPositionTSIdifference = abs(TruthTSI(2)-yPositionTSI)*1000; % y-position TSI [m]
% zPositionTSIdifference = abs(TruthTSI(3)-zPositionTSI)*1000; % z-position TSI [m]
% xVelocityTSIdifference = abs(TruthTSI(4)-xVelocityTSI)*1000; % x-velocity TSI [m/s]
% yVelocityTSIdifference = abs(TruthTSI(5)-yVelocityTSI)*1000; % y-velocity TSI [m/s]
% zVelocityTSIdifference = abs(TruthTSI(6)-zVelocityTSI)*1000; % z-velocity TSI [m/s]
% MassTSIdifference = abs(TruthTSI(7)-MassTSI); % Mass TSI [kg]
% 
% xPositionTSIdifferenceRKF = abs(VariableVectorTSI(:,21))*1000; % x-position TSI [m]
% yPositionTSIdifferenceRKF = abs(VariableVectorTSI(:,22))*1000; % y-position TSI [m]
% zPositionTSIdifferenceRKF = abs(VariableVectorTSI(:,23))*1000; % z-position TSI [m]
% xVelocityTSIdifferenceRKF = abs(VariableVectorTSI(:,24))*1000; % x-velocity TSI [m/s]
% yVelocityTSIdifferenceRKF = abs(VariableVectorTSI(:,25))*1000; % y-velocity TSI [m/s]
% zVelocityTSIdifferenceRKF = abs(VariableVectorTSI(:,26))*1000; % z-velocity TSI [m/s]
% MassTSIdifferenceRKF = abs(VariableVectorTSI(:,27)); % Mass TSI [kg]
% 



%% Plots

% x position
figure(1)
plot(timeTSIOrder20,xPositionTSIOrder20)
hold on
plot(timeTSIOrder29,xPositionTSIOrder29)
hold on
plot(timeTSIOrder30,xPositionTSIOrder30)

title(['Time vs ','x position all']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('x position [km]');

legend('Order 20','Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

xPositionDifference29 = xPositionTSIOrder20-xPositionTSIOrder29;
xPositionDifference30 = xPositionTSIOrder20-xPositionTSIOrder30;

figure(2)
plot(timeTSIOrder20,xPositionTSIOrder20-xPositionTSIOrder29)
hold on
plot(timeTSIOrder30,xPositionTSIOrder20-xPositionTSIOrder30)

title(['Time vs ','x position difference']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('x position difference [km]');

legend('Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

% figure(3)
% plot(timeTSIOrder20,xPositionTSIOrder20)
% hold on
% plot(timeTSIOrder29,xPositionTSIOrder29)
% 
% 
% title(['Time vs ','x position all faulty order']); % Give the figure a title
% xlabel('Time [s]'); % Label the different axes
% ylabel('x position [km]');
% 
% legend('Order 20','Order 29','Location','NorthEastOutside'); % Add a legend in the top right corner

% y position
figure(4)
plot(timeTSIOrder20,yPositionTSIOrder20)
hold on
plot(timeTSIOrder29,yPositionTSIOrder29)
hold on
plot(timeTSIOrder30,yPositionTSIOrder30)

title(['Time vs ','y position all']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('y position [km]');

legend('Order 20','Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

yPositionDifference29 = yPositionTSIOrder20-yPositionTSIOrder29;
yPositionDifference30 = yPositionTSIOrder20-yPositionTSIOrder30;

figure(5)
plot(timeTSIOrder20,yPositionTSIOrder20-yPositionTSIOrder29)
hold on
plot(timeTSIOrder30,yPositionTSIOrder20-yPositionTSIOrder30)

title(['Time vs ','y position difference']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('y position difference [km]');

legend('Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

% figure(6)
% plot(timeTSIOrder20,yPositionTSIOrder20)
% hold on
% plot(timeTSIOrder29,yPositionTSIOrder29)
% 
% 
% title(['Time vs ','y position all faulty order']); % Give the figure a title
% xlabel('Time [s]'); % Label the different axes
% ylabel('y position [km]');
% 
% legend('Order 20','Order 29','Location','NorthEastOutside'); % Add a legend in the top right corner



% z position
figure(7)
plot(timeTSIOrder20,zPositionTSIOrder20)
hold on
plot(timeTSIOrder29,zPositionTSIOrder29)
hold on
plot(timeTSIOrder30,zPositionTSIOrder30)

title(['Time vs ','z position all']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('z position [km]');

legend('Order 20','Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

zPositionDifference29 = zPositionTSIOrder20-zPositionTSIOrder29;
zPositionDifference30 = zPositionTSIOrder20-zPositionTSIOrder30;

figure(8)
plot(timeTSIOrder20,zPositionTSIOrder20-zPositionTSIOrder29)
hold on
plot(timeTSIOrder30,zPositionTSIOrder20-zPositionTSIOrder30)

title(['Time vs ','z position difference']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('z position difference [km]');

legend('Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

% figure(9)
% plot(timeTSIOrder20,zPositionTSIOrder20)
% hold on
% plot(timeTSIOrder29,zPositionTSIOrder29)
% 
% 
% title(['Time vs ','z position all faulty order']); % Give the figure a title
% xlabel('Time [s]'); % Label the different axes
% ylabel('z position [km]');
% 
% legend('Order 20','Order 29','Location','NorthEastOutside'); % Add a legend in the top right corner

% x velocity
figure(10)
plot(timeTSIOrder20,xVelocityTSIOrder20)
hold on
plot(timeTSIOrder29,xVelocityTSIOrder29)
hold on
plot(timeTSIOrder30,xVelocityTSIOrder30)

title(['Time vs ','x velocity all']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('x velocity [km]');

legend('Order 20','Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

xVelocityDifference29 = xVelocityTSIOrder20-xVelocityTSIOrder29;
xVelocityDifference30 = xVelocityTSIOrder20-xVelocityTSIOrder30;

figure(11)
plot(timeTSIOrder20,xVelocityTSIOrder20-xVelocityTSIOrder29)
hold on
plot(timeTSIOrder30,xVelocityTSIOrder20-xVelocityTSIOrder30)

title(['Time vs ','x velocity difference']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('x velocity difference [km]');

legend('Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

% figure(12)
% plot(timeTSIOrder20,xVelocityTSIOrder20)
% hold on
% plot(timeTSIOrder29,xVelocityTSIOrder29)
% 
% 
% title(['Time vs ','x velocity all faulty order']); % Give the figure a title
% xlabel('Time [s]'); % Label the different axes
% ylabel('x velocity [km]');
% 
% legend('Order 20','Order 29','Location','NorthEastOutside'); % Add a legend in the top right corner

% y velocity
figure(13)
plot(timeTSIOrder20,yVelocityTSIOrder20)
hold on
plot(timeTSIOrder29,yVelocityTSIOrder29)
hold on
plot(timeTSIOrder30,yVelocityTSIOrder30)

title(['Time vs ','y velocity all']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('y velocity [km]');

legend('Order 20','Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

yVelocityDifference29 = yVelocityTSIOrder20-yVelocityTSIOrder29;
yVelocityDifference30 = yVelocityTSIOrder20-yVelocityTSIOrder30;

figure(14)
plot(timeTSIOrder20,yVelocityTSIOrder20-yVelocityTSIOrder29)
hold on
plot(timeTSIOrder30,yVelocityTSIOrder20-yVelocityTSIOrder30)

title(['Time vs ','y velocity difference']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('y velocity difference [km]');

legend('Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

% figure(15)
% plot(timeTSIOrder20,yVelocityTSIOrder20)
% hold on
% plot(timeTSIOrder29,yVelocityTSIOrder29)
% 
% 
% title(['Time vs ','y velocity all faulty order']); % Give the figure a title
% xlabel('Time [s]'); % Label the different axes
% ylabel('y velocity [km]');
% 
% legend('Order 20','Order 29','Location','NorthEastOutside'); % Add a legend in the top right corner

% z velocity
figure(16)
plot(timeTSIOrder20,zVelocityTSIOrder20)
hold on
plot(timeTSIOrder29,zVelocityTSIOrder29)
hold on
plot(timeTSIOrder30,zVelocityTSIOrder30)

title(['Time vs ','z velocity all']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('z velocity [km]');

legend('Order 20','Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

zVelocityDifference29 = zVelocityTSIOrder20-zVelocityTSIOrder29;
zVelocityDifference30 = zVelocityTSIOrder20-zVelocityTSIOrder30;

figure(17)
plot(timeTSIOrder20,zVelocityTSIOrder20-zVelocityTSIOrder29)
hold on
plot(timeTSIOrder30,zVelocityTSIOrder20-zVelocityTSIOrder30)

title(['Time vs ','z velocity difference']); % Give the figure a title
xlabel('Time [s]'); % Label the different axes
ylabel('z velocity difference [km]');

legend('Order 29','Order 30','Location','NorthEastOutside'); % Add a legend in the top right corner

% figure(18)
% plot(timeTSIOrder20,zVelocityTSIOrder20)
% hold on
% plot(timeTSIOrder29,zVelocityTSIOrder29)
% 
% 
% title(['Time vs ','z velocity all faulty order']); % Give the figure a title
% xlabel('Time [s]'); % Label the different axes
% ylabel('z velocity [km]');
% 
% legend('Order 20','Order 29','Location','NorthEastOutside'); % Add a legend in the top right corner

toc




