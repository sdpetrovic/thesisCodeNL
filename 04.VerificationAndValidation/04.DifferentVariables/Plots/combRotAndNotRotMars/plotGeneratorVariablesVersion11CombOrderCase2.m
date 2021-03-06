% This script lets me plot the 3-D graphs of the generated trajectories
% and plot the sensitivity analysis plots.
% For this the Cartesian generated data is used (also in the spherical
% computation case) and the latests files for the different variables.
%
% Stacha Petrovic 15-12-2016
% 
% version 11
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

pathToValidationFolder = '/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/04.VerificationAndValidation/04.DifferentVariables/';

currentVariableFolder = 'Order';
currentVariableUnit = '[-]';
currentCaseFolder = 'Case2';
filePref = 0; % Choose a different file from the last one. So if this number is 1, you choose the penultimate file etc.

% Automatic latest file name (Source:
% https://www.youtube.com/watch?v=D8UkAOhsyDI)

% Write the paths to the folders
pathTSIrot = fullfile(pathToValidationFolder,'rotatingMars/TSI',currentVariableFolder,currentCaseFolder);
pathTSInotRot = fullfile(pathToValidationFolder,'nonRotatingMars/TSI',currentVariableFolder,currentCaseFolder);

pathRKF = fullfile(pathToValidationFolder,'RKF',currentVariableFolder,currentCaseFolder);

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
cpuTimeTSI_rot = VariableVectorTSIrot(:,3); % The CPU time estimate for TSI [sec]
wallTimeTSI_rot = VariableVectorTSIrot(:,4); % The wall time estimate for TSI [sec]
functionEvaluationsTSI_rot = VariableVectorTSIrot(:,5); % The TSI function evaluations

currentVariable_notRot = VariableVectorTSInotRot(:,1); % Current variable values
cpuTimeTSI_notRot = VariableVectorTSInotRot(:,3); % The CPU time estimate for TSI [sec]
wallTimeTSI_notRot = VariableVectorTSInotRot(:,4); % The wall time estimate for TSI [sec]
functionEvaluationsTSI_notRot = VariableVectorTSInotRot(:,5); % The TSI function evaluations


% End states
xPositionTSI_rot = VariableVectorTSIrot(:,7); % x-position TSI [km]
yPositionTSI_rot = VariableVectorTSIrot(:,8); % y-position TSI [km]
zPositionTSI_rot = VariableVectorTSIrot(:,9); % z-position TSI [km]
xVelocityTSI_rot = VariableVectorTSIrot(:,10); % x-velocity TSI [km/s]
yVelocityTSI_rot = VariableVectorTSIrot(:,11); % y-velocity TSI [km/s]
zVelocityTSI_rot = VariableVectorTSIrot(:,12); % z-velocity TSI [km/s]
MassTSI_rot = VariableVectorTSIrot(:,13); % Mass TSI [kg]

% Required propellant mass to circularize and reach desired inclination
propMassTSI_rot = VariableVectorTSIrot(:,28); % TSI required propellant mass [kg]


% End states
xPositionTSI_notRot = VariableVectorTSInotRot(:,7); % x-position TSI [km]
yPositionTSI_notRot = VariableVectorTSInotRot(:,8); % y-position TSI [km]
zPositionTSI_notRot = VariableVectorTSInotRot(:,9); % z-position TSI [km]
xVelocityTSI_notRot = VariableVectorTSInotRot(:,10); % x-velocity TSI [km/s]
yVelocityTSI_notRot = VariableVectorTSInotRot(:,11); % y-velocity TSI [km/s]
zVelocityTSI_notRot = VariableVectorTSInotRot(:,12); % z-velocity TSI [km/s]
MassTSI_notRot = VariableVectorTSInotRot(:,13); % Mass TSI [kg]

% Required propellant mass to circularize and reach desired inclination
propMassTSI_notRot = VariableVectorTSInotRot(:,28); % TSI required propellant mass [kg]


% Accuracy of the results and speed of convergence
% TruthTSIrot = VariableVectorTSIrot(1,7:13);
TruthTSIrot = VariableVectorTSIrot((find(VariableVectorTSIrot(:,1)==20)),7:13);
TruthTSInotRot = VariableVectorTSInotRot((find(VariableVectorTSInotRot(:,1)==20)),7:13);
index = (find(VariableVectorTSInotRot(:,1)==20));

% Difference in metres w.r.t. the order 20 end state

xPositionTSIdifference_rot = abs(TruthTSIrot(1)-xPositionTSI_rot)*1000; % x-position TSI [m]
xPositionTSIdifference_rot(index) = 10^(-15);
yPositionTSIdifference_rot = abs(TruthTSIrot(2)-yPositionTSI_rot)*1000; % y-position TSI [m]
yPositionTSIdifference_rot(index) = 10^(-15);
zPositionTSIdifference_rot = abs(TruthTSIrot(3)-zPositionTSI_rot)*1000; % z-position TSI [m]
zPositionTSIdifference_rot(index) = 10^(-15);
xVelocityTSIdifference_rot = abs(TruthTSIrot(4)-xVelocityTSI_rot)*1000; % x-velocity TSI [m/s]
xVelocityTSIdifference_rot(index) = 10^(-15);
yVelocityTSIdifference_rot = abs(TruthTSIrot(5)-yVelocityTSI_rot)*1000; % y-velocity TSI [m/s]
yVelocityTSIdifference_rot(index) = 10^(-15);
zVelocityTSIdifference_rot = abs(TruthTSIrot(6)-zVelocityTSI_rot)*1000; % z-velocity TSI [m/s]
zVelocityTSIdifference_rot(index) = 10^(-15);
MassTSIdifference_rot = abs(TruthTSIrot(7)-MassTSI_rot); % Mass TSI [kg]



xPositionTSIdifference_notRot = abs(TruthTSInotRot(1)-xPositionTSI_notRot)*1000; % x-position TSI [m]
xPositionTSIdifference_notRot(index) = 10^(-15);
yPositionTSIdifference_notRot = abs(TruthTSInotRot(2)-yPositionTSI_notRot)*1000; % y-position TSI [m]
yPositionTSIdifference_notRot(index) = 10^(-15);
zPositionTSIdifference_notRot = abs(TruthTSInotRot(3)-zPositionTSI_notRot)*1000; % z-position TSI [m]
zPositionTSIdifference_notRot(index) = 10^(-15);
xVelocityTSIdifference_notRot = abs(TruthTSInotRot(4)-xVelocityTSI_notRot)*1000; % x-velocity TSI [m/s]
xVelocityTSIdifference_notRot(index) = 10^(-15);
yVelocityTSIdifference_notRot = abs(TruthTSInotRot(5)-yVelocityTSI_notRot)*1000; % y-velocity TSI [m/s]
yVelocityTSIdifference_notRot(index) = 10^(-15);
zVelocityTSIdifference_notRot = abs(TruthTSInotRot(6)-zVelocityTSI_notRot)*1000; % z-velocity TSI [m/s]
zVelocityTSIdifference_notRot(index) = 10^(-15);
MassTSIdifference_notRot = abs(TruthTSInotRot(7)-MassTSI_notRot); % Mass TSI [kg]

xPositionTSIdifferenceRKF_rot = abs(VariableVectorTSIrot(:,21))*1000; % x-position TSI [m]
yPositionTSIdifferenceRKF_rot = abs(VariableVectorTSIrot(:,22))*1000; % y-position TSI [m]
zPositionTSIdifferenceRKF_rot = abs(VariableVectorTSIrot(:,23))*1000; % z-position TSI [m]
xVelocityTSIdifferenceRKF_rot = abs(VariableVectorTSIrot(:,24))*1000; % x-velocity TSI [m/s]
yVelocityTSIdifferenceRKF_rot = abs(VariableVectorTSIrot(:,25))*1000; % y-velocity TSI [m/s]
zVelocityTSIdifferenceRKF_rot = abs(VariableVectorTSIrot(:,26))*1000; % z-velocity TSI [m/s]
MassTSIdifferenceRKF_rot = abs(VariableVectorTSIrot(:,27)); % Mass TSI [kg]

xPositionTSIdifferenceRKF_notRot = abs(VariableVectorTSInotRot(:,21))*1000; % x-position TSI [m]
yPositionTSIdifferenceRKF_notRot = abs(VariableVectorTSInotRot(:,22))*1000; % y-position TSI [m]
zPositionTSIdifferenceRKF_notRot = abs(VariableVectorTSInotRot(:,23))*1000; % z-position TSI [m]
xVelocityTSIdifferenceRKF_notRot = abs(VariableVectorTSInotRot(:,24))*1000; % x-velocity TSI [m/s]
yVelocityTSIdifferenceRKF_notRot = abs(VariableVectorTSInotRot(:,25))*1000; % y-velocity TSI [m/s]
zVelocityTSIdifferenceRKF_notRot = abs(VariableVectorTSInotRot(:,26))*1000; % z-velocity TSI [m/s]
MassTSIdifferenceRKF_notRot = abs(VariableVectorTSInotRot(:,27)); % Mass TSI [kg]


% Convergence through difference between each consecutive value

for k = 1:1:length(xPositionTSI_rot)
   
    if k ~= length(xPositionTSI_rot)
    xPositionTSIdifferenceCons_rot(k) = abs(xPositionTSI_rot(k+1)-xPositionTSI_rot(k))*1000; % x-position TSI [m]
yPositionTSIdifferenceCons_rot(k) = abs(yPositionTSI_rot(k+1)-yPositionTSI_rot(k))*1000; % y-position TSI [m]
zPositionTSIdifferenceCons_rot(k) = abs(zPositionTSI_rot(k+1)-zPositionTSI_rot(k))*1000; % z-position TSI [m]
xVelocityTSIdifferenceCons_rot(k) = abs(xVelocityTSI_rot(k+1)-xVelocityTSI_rot(k))*1000; % x-velocity TSI [m/s]
yVelocityTSIdifferenceCons_rot(k) = abs(yVelocityTSI_rot(k+1)-yVelocityTSI_rot(k))*1000; % y-velocity TSI [m/s]
zVelocityTSIdifferenceCons_rot(k) = abs(zVelocityTSI_rot(k+1)-zVelocityTSI_rot(k))*1000; % z-velocity TSI [m/s]
MassTSIdifferenceCons_rot(k) = abs(MassTSI_rot(k+1)-MassTSI_rot(k)); % Mass TSI [kg]
    end

if k == length(xPositionTSI_rot)
   
        xPositionTSIdifferenceCons_rot(k) = 10^(-15); % x-position TSI [m]
yPositionTSIdifferenceCons_rot(k) = 10^(-15); % y-position TSI [m]
zPositionTSIdifferenceCons_rot(k) = 10^(-15); % z-position TSI [m]
xVelocityTSIdifferenceCons_rot(k) = 10^(-15); % x-velocity TSI [m/s]
yVelocityTSIdifferenceCons_rot(k) = 10^(-15); % y-velocity TSI [m/s]
zVelocityTSIdifferenceCons_rot(k) = 10^(-15); % z-velocity TSI [m/s]
MassTSIdifferenceCons_rot(k) = 10^(-15); % Mass TSI [kg]
    
end

end



for k = 1:1:length(xPositionTSI_notRot)
   
    if k ~= length(xPositionTSI_notRot)
    xPositionTSIdifferenceCons_notRot(k) = abs(xPositionTSI_notRot(k+1)-xPositionTSI_notRot(k))*1000; % x-position TSI [m]
yPositionTSIdifferenceCons_notRot(k) = abs(yPositionTSI_notRot(k+1)-yPositionTSI_notRot(k))*1000; % y-position TSI [m]
zPositionTSIdifferenceCons_notRot(k) = abs(zPositionTSI_notRot(k+1)-zPositionTSI_notRot(k))*1000; % z-position TSI [m]
xVelocityTSIdifferenceCons_notRot(k) = abs(xVelocityTSI_notRot(k+1)-xVelocityTSI_notRot(k))*1000; % x-velocity TSI [m/s]
yVelocityTSIdifferenceCons_notRot(k) = abs(yVelocityTSI_notRot(k+1)-yVelocityTSI_notRot(k))*1000; % y-velocity TSI [m/s]
zVelocityTSIdifferenceCons_notRot(k) = abs(zVelocityTSI_notRot(k+1)-zVelocityTSI_notRot(k))*1000; % z-velocity TSI [m/s]
MassTSIdifferenceCons_notRot(k) = abs(MassTSI_notRot(k+1)-MassTSI_notRot(k)); % Mass TSI [kg]
    end

if k == length(xPositionTSI_notRot)
   
        xPositionTSIdifferenceCons_notRot(k) = 10^(-15); % x-position TSI [m]
yPositionTSIdifferenceCons_notRot(k) = 10^(-15); % y-position TSI [m]
zPositionTSIdifferenceCons_notRot(k) = 10^(-15); % z-position TSI [m]
xVelocityTSIdifferenceCons_notRot(k) = 10^(-15); % x-velocity TSI [m/s]
yVelocityTSIdifferenceCons_notRot(k) = 10^(-15); % y-velocity TSI [m/s]
zVelocityTSIdifferenceCons_notRot(k) = 10^(-15); % z-velocity TSI [m/s]
MassTSIdifferenceCons_notRot(k) = 10^(-15); % Mass TSI [kg]
    
end

end


%% Plots

%% Evaluations per currentVariable and end state

figure(1) % CPU time
scatter(currentVariable_rot,cpuTimeTSI_rot);
hold on
scatter(currentVariable_notRot,cpuTimeTSI_notRot);


%title(['CPU time TSI vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('CPU time TSI [sec]');

axis([min(currentVariable_rot) max(currentVariable_rot) 0.0 0.05]) % First zoom
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom

% axis([currentVariable(1) currentVariable(end) 0.05 0.16]) % First zoom
% case 2
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom
% case 2


% legend('CPU time','Location','NorthEastOutside'); % Add a legend in the top right corner
legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(2) % Wall time
plot(currentVariable_rot,wallTimeTSI_rot);
hold on
plot(currentVariable_notRot,wallTimeTSI_notRot);


%title(['Wall time TSI vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Wall time TSI [sec]');

axis([min(currentVariable_rot) max(currentVariable_rot) 0.0 0.05]) % First zoom
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom

% axis([currentVariable(1) currentVariable(end) 0.05 0.16]) % First zoom
% case 2
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom
% case 2

% legend('Wall time','Location','NorthEastOutside'); % Add a legend in the top right corner
legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(3) % Position
plot(currentVariable_rot,xPositionTSI_rot)
hold on
plot(currentVariable_notRot,xPositionTSI_notRot)


%title(['x-Position vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('x-Position [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean([mean(xPositionTSI_rot),mean(xPositionTSI_notRot)]); % Mean of the vector
STD = std(xPositionTSI_rot);      % Standard deviation of the vector
dummyVector = xPositionTSI_rot;    % Create a replacement vector that can be adjusted

j = 1; % Used to determine the position of the value that has to be removed
for i = 1:1:length(dummyVector) % Replace every outlier that is beyond the standard deviation of the actual vector with the mean
    
        if dummyVector(i) > (Mean+STD)   % This is done to get an estimate of a proper axis boundary
        pos(j) = i;
        j = j+1;
    elseif dummyVector(i) < (Mean - STD)
        pos(j) = i;   
        j = j+1;
        end
                    
end

dummyVector(pos) = [];  % Delete the outliers

Meanadj = mean(dummyVector); % Determine the new mean (or adjusted mean, adj)
STDadj = std(dummyVector); % Determine the new std (or adjusted, adj)

% axis([currentVariable_rot(1) currentVariable_rot(end) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom



% legend('x-Position','y-Position','z-Position','Location','NorthEastOutside'); % Add a legend in the top right corner
legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(4)
plot(currentVariable_rot,yPositionTSI_rot)
hold on
plot(currentVariable_notRot,yPositionTSI_notRot)


%title(['y-Position vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('y-Position [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean([mean(yPositionTSI_rot),mean(yPositionTSI_notRot)]); % Mean of the vector
STD = std(yPositionTSI_rot);      % Standard deviation of the vector
dummyVector = yPositionTSI_rot;    % Create a replacement vector that can be adjusted

j = 1; % Used to determine the position of the value that has to be removed
for i = 1:1:length(dummyVector) % Replace every outlier that is beyond the standard deviation of the actual vector with the mean
    
        if dummyVector(i) > (Mean+STD)   % This is done to get an estimate of a proper axis boundary
        pos(j) = i;
        j = j+1;
    elseif dummyVector(i) < (Mean - STD)
        pos(j) = i;   
        j = j+1;
        end
                    
end

dummyVector(pos) = [];  % Delete the outliers

Meanadj = mean(dummyVector); % Determine the new mean (or adjusted mean, adj)
STDadj = std(dummyVector); % Determine the new std (or adjusted, adj)

% axis([currentVariable_rot(1) currentVariable_rot(end) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Position','y-Position','z-Position','Location','NorthEastOutside'); % Add a legend in the top right corner
legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(5)
plot(currentVariable_rot,zPositionTSI_rot)
hold on
plot(currentVariable_notRot,zPositionTSI_notRot)


%title(['z-Position vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('z-Position [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean([mean(zPositionTSI_rot),mean(zPositionTSI_notRot)]); % Mean of the vector
STD = std(zPositionTSI_rot);      % Standard deviation of the vector
dummyVector = zPositionTSI_rot;    % Create a replacement vector that can be adjusted

j = 1; % Used to determine the position of the value that has to be removed
for i = 1:1:length(dummyVector) % Replace every outlier that is beyond the standard deviation of the actual vector with the mean
    
        if dummyVector(i) > (Mean+STD)   % This is done to get an estimate of a proper axis boundary
        pos(j) = i;
        j = j+1;
    elseif dummyVector(i) < (Mean - STD)
        pos(j) = i;   
        j = j+1;
        end
                    
end

dummyVector(pos) = [];  % Delete the outliers

Meanadj = mean(dummyVector); % Determine the new mean (or adjusted mean, adj)
STDadj = std(dummyVector); % Determine the new std (or adjusted, adj)

% axis([currentVariable_rot(1) currentVariable_rot(end) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Position','y-Position','z-Position','Location','NorthEastOutside'); % Add a legend in the top right corner
legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(6) % Velocity
plot(currentVariable_rot,xVelocityTSI_rot)
hold on
plot(currentVariable_notRot,xVelocityTSI_notRot)

% plot(currentVariable_rot,yVelocityTSI)
% hold on
% plot(currentVariable_rot,zVelocityTSI)

%title(['x-Velocity vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('x-Velocity [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean([mean(xVelocityTSI_rot),mean(xVelocityTSI_notRot)]); % Mean of the vector
STD = std(xVelocityTSI_rot);      % Standard deviation of the vector
dummyVector = xVelocityTSI_rot;    % Create a replacement vector that can be adjusted

j = 1; % Used to determine the position of the value that has to be removed
for i = 1:1:length(dummyVector) % Replace every outlier that is beyond the standard deviation of the actual vector with the mean
    
        if dummyVector(i) > (Mean+STD)   % This is done to get an estimate of a proper axis boundary
        pos(j) = i;
        j = j+1;
    elseif dummyVector(i) < (Mean - STD)
        pos(j) = i;   
        j = j+1;
        end
                    
end

dummyVector(pos) = [];  % Delete the outliers

Meanadj = mean(dummyVector); % Determine the new mean (or adjusted mean, adj)
STDadj = std(dummyVector); % Determine the new std (or adjusted, adj)

% axis([currentVariable_rot(1) currentVariable_rot(end) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Velocity','y-Velocity','z-Velocity','Location','NorthEastOutside'); % Add a legend in the top right corner
legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(7)
plot(currentVariable_rot,yVelocityTSI_rot)
hold on
plot(currentVariable_notRot,yVelocityTSI_notRot)


%title(['y-Velocity vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('y-Velocity [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean([mean(yVelocityTSI_rot),mean(yVelocityTSI_notRot)]); % Mean of the vector
STD = std(yVelocityTSI_rot);      % Standard deviation of the vector
dummyVector = yVelocityTSI_rot;    % Create a replacement vector that can be adjusted

j = 1; % Used to determine the position of the value that has to be removed
for i = 1:1:length(dummyVector) % Replace every outlier that is beyond the standard deviation of the actual vector with the mean
    
        if dummyVector(i) > (Mean+STD)   % This is done to get an estimate of a proper axis boundary
        pos(j) = i;
        j = j+1;
    elseif dummyVector(i) < (Mean - STD)
        pos(j) = i;   
        j = j+1;
        end
                    
end

dummyVector(pos) = [];  % Delete the outliers

Meanadj = mean(dummyVector); % Determine the new mean (or adjusted mean, adj)
STDadj = std(dummyVector); % Determine the new std (or adjusted, adj)

% axis([currentVariable_rot(1) currentVariable_rot(end) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Velocity','y-Velocity','z-Velocity','Location','NorthEastOutside'); % Add a legend in the top right corner
legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(8)
plot(currentVariable_rot,zVelocityTSI_rot)
hold on
plot(currentVariable_notRot,zVelocityTSI_notRot)


%title(['z-Velocity vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('z-Velocity [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean([mean(zVelocityTSI_rot),mean(zVelocityTSI_notRot)]); % Mean of the vector
STD = std(zVelocityTSI_rot);      % Standard deviation of the vector
dummyVector = zVelocityTSI_rot;    % Create a replacement vector that can be adjusted

j = 1; % Used to determine the position of the value that has to be removed
for i = 1:1:length(dummyVector) % Replace every outlier that is beyond the standard deviation of the actual vector with the mean
    
        if dummyVector(i) > (Mean+STD)   % This is done to get an estimate of a proper axis boundary
        pos(j) = i;
        j = j+1;
    elseif dummyVector(i) < (Mean - STD)
        pos(j) = i;   
        j = j+1;
        end
                    
end

dummyVector(pos) = [];  % Delete the outliers

Meanadj = mean(dummyVector); % Determine the new mean (or adjusted mean, adj)
STDadj = std(dummyVector); % Determine the new std (or adjusted, adj)

% axis([currentVariable_rot(1) currentVariable_rot(end) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Velocity','y-Velocity','z-Velocity','Location','NorthEastOutside'); % Add a legend in the top right corner
legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14


figure(9)
semilogy(currentVariable_rot,functionEvaluationsTSI_rot)
hold on
semilogy(currentVariable_notRot,functionEvaluationsTSI_notRot)


%title(['Function evaluations vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Function evaluations [-]');

% Determine the mean and the standard deviation of the vector
Mean = mean(functionEvaluationsTSI_rot); % Mean of the vector
STD = std(functionEvaluationsTSI_rot);      % Standard deviation of the vector
dummyVector = functionEvaluationsTSI_rot;    % Create a replacement vector that can be adjusted

j = 1; % Used to determine the position of the value that has to be removed
for i = 1:1:length(dummyVector) % Replace every outlier that is beyond the standard deviation of the actual vector with the mean
    
        if dummyVector(i) > (Mean+STD)   % This is done to get an estimate of a proper axis boundary
        pos(j) = i;
        j = j+1;
    elseif dummyVector(i) < (Mean - STD)
        pos(j) = i;   
        j = j+1;
        end
                    
end

dummyVector(pos) = [];  % Delete the outliers

Meanadj = mean(dummyVector); % Determine the new mean (or adjusted mean, adj)
STDadj = std(dummyVector); % Determine the new std (or adjusted, adj)

% axis([currentVariable(1) currentVariable(end) 0 Meanadj+3*STDadj]) % Zoom
% axis([currentVariable(1) currentVariable(end) 0 100]) % Zoom

legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(10)
% semilogx(functionEvaluationsTSI,cpuTimeTSI)
scatter(functionEvaluationsTSI_rot,cpuTimeTSI_rot)
hold on
scatter(functionEvaluationsTSI_notRot,cpuTimeTSI_notRot)

set(gca,'xscale','log')



%title(['Function evaluations vs CPU time TSI']); % Give the figure a title
xlabel('Function evaluations [-]');
ylabel('CPU time TSI [sec]'); % Label the different axes


% axis([min(functionEvaluationsTSI) max(functionEvaluationsTSI) 0.05 0.09]) % Zoom
axis([min(functionEvaluationsTSI_rot) max(functionEvaluationsTSI_rot) 0.0 0.04]) % Zoom

legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
set(gca,'FontSize',14);  % Make sure that the graph font size is 14

%% Speed of convergence and accuracy

figure(11)
semilogy(currentVariable_rot,xPositionTSIdifference_rot)
hold on
semilogy(currentVariable_rot,yPositionTSIdifference_rot)
hold on
semilogy(currentVariable_rot,zPositionTSIdifference_rot)
hold on
semilogy(currentVariable_rot,xVelocityTSIdifference_rot)
hold on
semilogy(currentVariable_rot,yVelocityTSIdifference_rot)
hold on
semilogy(currentVariable_rot,zVelocityTSIdifference_rot)
hold on 
% scatter(currentVariable_notRot,xPositionTSIdifference_notRot)
% hold on
% semilogy(currentVariable_rot,MassTSIdifference_rot)


%title(['Absolute differences rotating Mars w.r.t. the nominal case for every ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Absolute difference [m, m/s]');

legend('x-Position rot','y-Position rot','z-Position rot','x-Velocity rot','y-Velocity rot','z-Velocity rot','x-Position notRot','y-Position notRot','z-Position notRot','x-Velocity notRot','y-Velocity notRot','z-Velocity notRot','Location','NorthEastOutside'); % Add a legend in the top right corner

set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(12)
semilogy(currentVariable_rot,xPositionTSIdifferenceRKF_rot)
hold on
semilogy(currentVariable_rot,yPositionTSIdifferenceRKF_rot)
hold on
semilogy(currentVariable_rot,zPositionTSIdifferenceRKF_rot)
hold on
semilogy(currentVariable_rot,xVelocityTSIdifferenceRKF_rot)
hold on
semilogy(currentVariable_rot,yVelocityTSIdifferenceRKF_rot)
hold on
semilogy(currentVariable_rot,zVelocityTSIdifferenceRKF_rot)
hold on
% semilogy(currentVariable_rot,MassTSIdifferenceRKF_rot)


%title(['Absolute differences rotating Mars w.r.t. the RKF case for every ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Absolute difference [m, m/s]');

legend('x-Position rot','y-Position rot','z-Position rot','x-Velocity rot','y-Velocity rot','z-Velocity rot','x-Position notRot','y-Position notRot','z-Position notRot','x-Velocity notRot','y-Velocity notRot','z-Velocity notRot','Location','NorthEastOutside'); % Add a legend in the top right corner

set(gca,'FontSize',14);  % Make sure that the graph font size is 14



figure(13)
semilogy(currentVariable_rot,xPositionTSIdifferenceCons_rot)
hold on
semilogy(currentVariable_rot,yPositionTSIdifferenceCons_rot)
hold on
semilogy(currentVariable_rot,zPositionTSIdifferenceCons_rot)
hold on
semilogy(currentVariable_rot,xVelocityTSIdifferenceCons_rot)
hold on
semilogy(currentVariable_rot,yVelocityTSIdifferenceCons_rot)
hold on
semilogy(currentVariable_rot,zVelocityTSIdifferenceCons_rot)
hold on
% semilogy(currentVariable_rot,MassTSIdifferenceCons_rot)

xlim([0 30]);

%title(['Consecutive differences rotating Mars between ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Consecutive differences [m, m/s]');

legend('x-Position rot','y-Position rot','z-Position rot','x-Velocity rot','y-Velocity rot','z-Velocity rot','x-Position notRot','y-Position notRot','z-Position notRot','x-Velocity notRot','y-Velocity notRot','z-Velocity notRot','Location','NorthEastOutside'); % Add a legend in the top right corner

set(gca,'FontSize',14);  % Make sure that the graph font size is 14


%%%%%%%%%%%%% Non rotating

figure(11)
scatter(currentVariable_notRot,xPositionTSIdifference_notRot)
hold on
scatter(currentVariable_notRot,yPositionTSIdifference_notRot)
hold on
scatter(currentVariable_notRot,zPositionTSIdifference_notRot)
hold on
scatter(currentVariable_notRot,xVelocityTSIdifference_notRot)
hold on
scatter(currentVariable_notRot,yVelocityTSIdifference_notRot)
hold on
scatter(currentVariable_notRot,zVelocityTSIdifference_notRot)
% hold on 
% scatter(currentVariable_notnotRot,xPositionTSIdifference_notnotRot)
% hold on
% scatter(currentVariable_notRot,MassTSIdifference_notRot)


%title(['Absolute differences non-rotating Mars w.r.t. the nominal case for every ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Absolute difference [m, m/s]');

legend('x-Position rot','y-Position rot','z-Position rot','x-Velocity rot','y-Velocity rot','z-Velocity rot','x-Position notRot','y-Position notRot','z-Position notRot','x-Velocity notRot','y-Velocity notRot','z-Velocity notRot','Location','NorthEastOutside'); % Add a legend in the top right corner

set(gca,'FontSize',14);  % Make sure that the graph font size is 14

figure(12)
scatter(currentVariable_notRot,xPositionTSIdifferenceRKF_notRot)
hold on
scatter(currentVariable_notRot,yPositionTSIdifferenceRKF_notRot)
hold on
scatter(currentVariable_notRot,zPositionTSIdifferenceRKF_notRot)
hold on
scatter(currentVariable_notRot,xVelocityTSIdifferenceRKF_notRot)
hold on
scatter(currentVariable_notRot,yVelocityTSIdifferenceRKF_notRot)
hold on
scatter(currentVariable_notRot,zVelocityTSIdifferenceRKF_notRot)
% hold on
% scatter(currentVariable_notRot,MassTSIdifferenceRKF_notRot)


%title(['Absolute differences non-rotating Mars w.r.t. the RKF case for every ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Absolute difference [m, m/s]');

legend('x-Position rot','y-Position rot','z-Position rot','x-Velocity rot','y-Velocity rot','z-Velocity rot','x-Position notRot','y-Position notRot','z-Position notRot','x-Velocity notRot','y-Velocity notRot','z-Velocity notRot','Location','NorthEastOutside'); % Add a legend in the top right corner

set(gca,'FontSize',14);  % Make sure that the graph font size is 14



figure(13)
scatter(currentVariable_notRot,xPositionTSIdifferenceCons_notRot)
hold on
scatter(currentVariable_notRot,yPositionTSIdifferenceCons_notRot)
hold on
scatter(currentVariable_notRot,zPositionTSIdifferenceCons_notRot)
hold on
scatter(currentVariable_notRot,xVelocityTSIdifferenceCons_notRot)
hold on
scatter(currentVariable_notRot,yVelocityTSIdifferenceCons_notRot)
hold on
scatter(currentVariable_notRot,zVelocityTSIdifferenceCons_notRot)
% hold on
% scatter(currentVariable_rot,MassTSIdifferenceCons_rot)


%title(['Consecutive differences non-rotating Mars between ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Consecutive differences [m, m/s]');

xlim([0 30]);

legend('x-Position rot','y-Position rot','z-Position rot','x-Velocity rot','y-Velocity rot','z-Velocity rot','x-Position notRot','y-Position notRot','z-Position notRot','x-Velocity notRot','y-Velocity notRot','z-Velocity notRot','Location','NorthEastOutside'); % Add a legend in the top right corner

set(gca,'FontSize',14);  % Make sure that the graph font size is 14

% figure(17) % CPU time
% scatter(currentVariable_rot,cpuTimeTSI_rot);
% hold on
% 
% %%title(['CPU time TSI vs ',currentVariableFolder]); % Give the figure a title
% xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
% ylabel('CPU time TSI [sec]');
% 
% axis([min(currentVariable_rot) max(currentVariable_rot) 0.0 0.05]) % First zoom
% 
% legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner
% set(gca,'FontSize',14);  % Make sure that the graph font size is 14
% 
% figure(17)
%  scatter(currentVariable_notRot,cpuTimeTSI_notRot);
%  legend('Rotating','Non rotating','Location','NorthEastOutside'); % Add a legend in the top right corner



%%%%% Combined graph data

figure(14) % CPU time
yyaxis left
plot(currentVariable_rot,cpuTimeTSI_rot);


xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('CPU time TSI [sec]');

axis([min(currentVariable_rot) max(currentVariable_rot) 0.0 0.05]) % No live data collection
% axis([min(currentVariable_rot)  max(currentVariable_rot) 0.08 0.16]) % Live data
% collection



set(gca,'FontSize',14);  % Make sure that the graph font size is 14


yyaxis right
semilogy(currentVariable_rot,functionEvaluationsTSI_rot)


ylabel('Function evaluations [-]');

toc




