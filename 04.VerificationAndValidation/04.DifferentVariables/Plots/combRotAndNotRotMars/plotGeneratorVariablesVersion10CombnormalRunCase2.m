% This script lets me plot the 3-D graphs of the generated trajectories
% and plot the sensitivity analysis plots.
% For this the Cartesian generated data is used (also in the spherical
% computation case) and the latests files for the different variables.
%
% Stacha Petrovic 15-12-2016
% 
% version 10
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

pathToValidationFolder = '/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/04.VerificationAndValidation/04.DifferentVariables';

currentVariableFolder = 'normalRun';

%     /* List of possible variable names
%      *
%      *  - Order
%      *  - ErrorTolerance
%      *  - launchAltitude
%      *  - launchLatitude
%      *  - launchLongitude
%      *  - Position         (Such as (0,0) (0,90) and (90,0) latitude and longitude)
%      *  - FPA
%      *  - headingAngle
%      *  - groundVelocity
%      *  - normalRun



currentVariableUnit = '[-]';
currentCaseFolder = 'Case2';
filePref = 1; % Choose a different file from the last one. So if this number is 1, you choose the penultimate file etc.

% Automatic latest file name (Source:
% https://www.youtube.com/watch?v=D8UkAOhsyDI)

% Write the paths to the folders
pathTSIrot = fullfile(pathToValidationFolder,'rotatingMars/TSI',currentVariableFolder,currentCaseFolder);
pathTSInotRot = fullfile(pathToValidationFolder,'nonRotatingMars/TSI',currentVariableFolder,currentCaseFolder);

% pathRKF = fullfile(pathToValidationFolder,'RKF',currentVariableFolder,currentCaseFolder);

% Find all files ending with .csv in that folder
dTSIrot = dir(fullfile(pathTSIrot,'*.csv'));
dTSInotRot = dir(fullfile(pathTSInotRot,'*.csv'));
% dRKF = dir(fullfile(pathRKF,'*.csv'));

% Find the latest file
datesTSIrot = [dTSIrot.datenum]; % Get all the dates
[~,newestIndexTSIrot] = max(datesTSIrot); % Get the index of the newest date
newestFileTSIrot = dTSIrot(newestIndexTSIrot-filePref); % Get the data corresponding to the newest file
newestFileNameTSIrot = newestFileTSIrot.name; % Get the name of the newest file

datesTSInotRot = [dTSInotRot.datenum]; % Get all the dates
[~,newestIndexTSInotRot] = max(datesTSInotRot); % Get the index of the newest date
newestFileTSInotRot = dTSInotRot(newestIndexTSInotRot-filePref); % Get the data corresponding to the newest file
newestFileNameTSInotRot = newestFileTSInotRot.name; % Get the name of the newest file

% datesRKF = [dRKF.datenum]; % Get all the dates
% [~,newestIndexRKF] = max(datesRKF); % Get the index of the newest date
% newestFileRKF = dRKF(newestIndexRKF-filePref); % Get the data corresponding to the newest file
% newestFileNameRKF = newestFileRKF.name; % Get the name of the newest file



% Specify the paths
VariablePathTSIrot = fullfile(pathTSIrot,newestFileNameTSIrot); % Create the path for the TSI data
VariablePathTSInotRot = fullfile(pathTSInotRot,newestFileNameTSInotRot); % Create the path for the TSI data

% VariablePathRKF = fullfile(pathRKF,newestFileNameRKF); % Create the path for the RKF data


VariableVectorTSIrot = csvread(VariablePathTSIrot); % Read the file
VariableVectorTSInotRot = csvread(VariablePathTSInotRot); % Read the file
% VariableVectorRKF = csvread(VariablePathRKF); % Read the file


% Data in the file:

%     /// Storing the output TSI
% 
%     outputMatrixTSI(run,0) = run+1;
% //    outputMatrixTSI(run,0) = initialLatitudeDeg;
% //    outputMatrixTSI(run,0) = FlightPathAngleDeg;
% //    outputMatrixTSI(run,0) = HeadingAngleDeg;
% //    outputMatrixTSI(run,0) = initialGroundVelocity;
% //    outputMatrixTSI(run,0) = maxOrder;
% //    outputMatrixTSI(run,1) = initialLongitudeDeg;
% 
%     outputMatrixTSI(run,2) = outputMatrix(3,0); // CPU time TSI
%     outputMatrixTSI(run,3) = outputMatrix(3,1); // Wall time TSI
%     outputMatrixTSI(run,4) = outputMatrix(3,4); // Number of evaluations TSI
% 
%     // State and time
%     outputMatrixTSI(run,5) = outputMatrix(0,0); // Time (End time)
%     outputMatrixTSI(run,6) = outputMatrix(0,1);
%     outputMatrixTSI(run,7) = outputMatrix(0,2);
%     outputMatrixTSI(run,8) = outputMatrix(0,3);
%     outputMatrixTSI(run,9) = outputMatrix(0,4);
%     outputMatrixTSI(run,10) = outputMatrix(0,5);
%     outputMatrixTSI(run,11) = outputMatrix(0,6);
%     outputMatrixTSI(run,12) = outputMatrix(0,7); // Mass
% 
%     // Fraction difference w.r.t. RKF
%     outputMatrixTSI(run,13) = outputMatrix(4,0);
%     outputMatrixTSI(run,14) = outputMatrix(4,1);
%     outputMatrixTSI(run,15) = outputMatrix(4,2);
%     outputMatrixTSI(run,16) = outputMatrix(4,3);
%     outputMatrixTSI(run,17) = outputMatrix(4,4);
%     outputMatrixTSI(run,18) = outputMatrix(4,5);
%     outputMatrixTSI(run,19) = outputMatrix(4,6);
% 
%     // Value difference w.r.t. RKF
%     outputMatrixTSI(run,20) = outputMatrix(5,0);
%     outputMatrixTSI(run,21) = outputMatrix(5,1);
%     outputMatrixTSI(run,22) = outputMatrix(5,2);
%     outputMatrixTSI(run,23) = outputMatrix(5,3);
%     outputMatrixTSI(run,24) = outputMatrix(5,4);
%     outputMatrixTSI(run,25) = outputMatrix(5,5);
%     outputMatrixTSI(run,26) = outputMatrix(5,6);
% 
%     // In case of circularisation:
%     outputMatrixTSI(run,27) = outputMatrix(2,4); // The propellant mass used for the circularisation burn of TSI
% 
%     /// Storing the output RKF
% 
%     outputMatrixRKF(run,0) = outputMatrixTSI(run,0);
%     outputMatrixRKF(run,1) = outputMatrixTSI(run,1);
% 
%     outputMatrixRKF(run,2) = outputMatrix(3,2); // CPU time RKF
%     outputMatrixRKF(run,3) = outputMatrix(3,3); // Wall time RKF
%     outputMatrixRKF(run,4) = outputMatrix(3,6); // Number of function evaluations RKF
% 
%     // State and time
%     outputMatrixRKF(run,5) = outputMatrix(1,0); // Time
%     outputMatrixRKF(run,6) = outputMatrix(1,1);
%     outputMatrixRKF(run,7) = outputMatrix(1,2);
%     outputMatrixRKF(run,8) = outputMatrix(1,3);
%     outputMatrixRKF(run,9) = outputMatrix(1,4);
%     outputMatrixRKF(run,10) = outputMatrix(1,5);
%     outputMatrixRKF(run,11) = outputMatrix(1,6);
%     outputMatrixRKF(run,12) = outputMatrix(1,7); // Mass
% 
%     // In case of circularisation:
%     outputMatrixRKF(run,13) = outputMatrix(2,5); // The propellant mass used for the circularisation burn of RKF




%% Generate the required vectors

% Desired data
currentVariable_rot = VariableVectorTSIrot(:,1); % Current variable values
Longitude_rot = VariableVectorTSIrot(:,2); % Initial longitude [deg]
cpuTimeTSIrot = VariableVectorTSIrot(:,3); % The CPU time estimate for TSI [sec]
wallTimeTSIrot = VariableVectorTSIrot(:,4); % The wall time estimate for TSI [sec]
functionEvaluationsTSIrot = VariableVectorTSIrot(:,5); % The TSI function evaluations

currentVariable_notRot = VariableVectorTSInotRot(:,1); % Current variable values
cpuTimeTSInotRot = VariableVectorTSInotRot(:,3); % The CPU time estimate for TSI [sec]
wallTimeTSInotRot = VariableVectorTSInotRot(:,4); % The wall time estimate for TSI [sec]
functionEvaluationsTSInotRot = VariableVectorTSInotRot(:,5); % The TSI function evaluations


% End states
xPositionTSIrot = VariableVectorTSIrot(:,7); % x-position TSI [km]
yPositionTSIrot = VariableVectorTSIrot(:,8); % y-position TSI [km]
zPositionTSIrot = VariableVectorTSIrot(:,9); % z-position TSI [km]
xVelocityTSIrot = VariableVectorTSIrot(:,10); % x-velocity TSI [km/s]
yVelocityTSIrot = VariableVectorTSIrot(:,11); % y-velocity TSI [km/s]
zVelocityTSIrot = VariableVectorTSIrot(:,12); % z-velocity TSI [km/s]
MassTSIrot = VariableVectorTSIrot(:,13); % Mass TSI [kg]

% Required propellant mass to circularize and reach desired inclination
propMassTSIrot = VariableVectorTSIrot(:,28); % TSI required propellant mass [kg]


% Accuracy of the results and speed of convergence
TruthTSIrot = VariableVectorTSIrot(1,7:13);
% TruthTSI = VariableVectorRKF((find(VariableVectorRKF(:,1)==20)),7:13);

% Difference in metres w.r.t. the order 20 end state

xPositionTSIdifference = abs(TruthTSIrot(1)-xPositionTSIrot)*1000; % x-position TSI [m]
yPositionTSIdifference = abs(TruthTSIrot(2)-yPositionTSIrot)*1000; % y-position TSI [m]
zPositionTSIdifference = abs(TruthTSIrot(3)-zPositionTSIrot)*1000; % z-position TSI [m]
xVelocityTSIdifference = abs(TruthTSIrot(4)-xVelocityTSIrot)*1000; % x-velocity TSI [m/s]
yVelocityTSIdifference = abs(TruthTSIrot(5)-yVelocityTSIrot)*1000; % y-velocity TSI [m/s]
zVelocityTSIdifference = abs(TruthTSIrot(6)-zVelocityTSIrot)*1000; % z-velocity TSI [m/s]
MassTSIdifference = abs(TruthTSIrot(7)-MassTSIrot); % Mass TSI [kg]

% xPositionTSIdifferenceRKF = abs(VariableVectorTSIrot(:,21))*1000; % x-position TSI [m]
% yPositionTSIdifferenceRKF = abs(VariableVectorTSIrot(:,22))*1000; % y-position TSI [m]
% zPositionTSIdifferenceRKF = abs(VariableVectorTSIrot(:,23))*1000; % z-position TSI [m]
% xVelocityTSIdifferenceRKF = abs(VariableVectorTSIrot(:,24))*1000; % x-velocity TSI [m/s]
% yVelocityTSIdifferenceRKF = abs(VariableVectorTSIrot(:,25))*1000; % y-velocity TSI [m/s]
% zVelocityTSIdifferenceRKF = abs(VariableVectorTSIrot(:,26))*1000; % z-velocity TSI [m/s]
% MassTSIdifferenceRKF = abs(VariableVectorTSIrot(:,27)); % Mass TSI [kg]

% Convergence through difference between each consecutive value

for k = 1:1:length(xPositionTSIrot)
   
    if k ~= length(xPositionTSIrot)
    xPositionTSIdifferenceCons(k) = abs(xPositionTSIrot(k+1)-xPositionTSIrot(k))*1000; % x-position TSI [m]
yPositionTSIdifferenceCons(k) = abs(yPositionTSIrot(k+1)-yPositionTSIrot(k))*1000; % y-position TSI [m]
zPositionTSIdifferenceCons(k) = abs(zPositionTSIrot(k+1)-zPositionTSIrot(k))*1000; % z-position TSI [m]
xVelocityTSIdifferenceCons(k) = abs(xVelocityTSIrot(k+1)-xVelocityTSIrot(k))*1000; % x-velocity TSI [m/s]
yVelocityTSIdifferenceCons(k) = abs(yVelocityTSIrot(k+1)-yVelocityTSIrot(k))*1000; % y-velocity TSI [m/s]
zVelocityTSIdifferenceCons(k) = abs(zVelocityTSIrot(k+1)-zVelocityTSIrot(k))*1000; % z-velocity TSI [m/s]
MassTSIdifferenceCons(k) = abs(MassTSIrot(k+1)-MassTSIrot(k)); % Mass TSI [kg]
    end

if k == length(xPositionTSIrot)
   
        xPositionTSIdifferenceCons(k) = 0; % x-position TSI [m]
yPositionTSIdifferenceCons(k) = 0; % y-position TSI [m]
zPositionTSIdifferenceCons(k) = 0; % z-position TSI [m]
xVelocityTSIdifferenceCons(k) = 0; % x-velocity TSI [m/s]
yVelocityTSIdifferenceCons(k) = 0; % y-velocity TSI [m/s]
zVelocityTSIdifferenceCons(k) = 0; % z-velocity TSI [m/s]
MassTSIdifferenceCons(k) = 0; % Mass TSI [kg]
    
end

end



%% Plots

%% Evaluations per currentVariable and end state

figure(1) % CPU time
% plot(currentVariable,cpuTimeTSI);
scatter(currentVariable_rot,cpuTimeTSIrot);
hold on
scatter(currentVariable_notRot,cpuTimeTSInotRot);

title(['CPU time TSI vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('CPU time TSI [sec]');

% axis([currentVariable(1) currentVariable(end) 0.05 0.1]) % First zoom
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom

% axis([currentVariable(1) currentVariable(end) 0.05 0.16]) % First zoom
% case 2
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom
% case 2


legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14



figure(2) % Wall time
plot(currentVariable_rot,wallTimeTSIrot);
hold on
plot(currentVariable_notRot,wallTimeTSInotRot);

title(['Wall time TSI vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Wall time TSI [sec]');

% axis([currentVariable_rot(1) currentVariable_rot(end) 0 max(max(wallTimeTSIrot),max(wallTimeTSInotRot))]) % First zoom
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom

% axis([currentVariable(1) currentVariable(end) 0.05 0.16]) % First zoom
% case 2
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom
% case 2

legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(10)
scatter(functionEvaluationsTSIrot,cpuTimeTSIrot)
hold on
myPoly = polyfit(functionEvaluationsTSIrot,cpuTimeTSIrot,3);
polyRange = min(functionEvaluationsTSIrot):1:max(functionEvaluationsTSIrot);
polyvalNumber = polyval(myPoly,polyRange);
hold on
scatter(functionEvaluationsTSInotRot,cpuTimeTSInotRot)

% plot(polyRange,polyvalNumber)

% plot(functionEvaluationsTSI,cpuTimeTSI)

title(['Function evaluations vs CPU time TSI']); % Give the figure a title
xlabel('Function evaluations [-]');
ylabel('CPU time TSI [sec]'); % Label the different axes

% xbounds = xlim();
% set(gca, 'xtick', xbounds(1):1:xbounds(2));



% axis([min(functionEvaluationsTSI)-5 100 min(cpuTimeTSI) max(cpuTimeTSI)]) % Zoom

legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

toc




