% This script lets me plot the 3-D graphs of the generated trajectories
% and plot the sensitivity analysis plots.
% For this the Cartesian generated data is used (also in the spherical
% computation case) and the latests files for the different variables.
%
% Stacha Petrovic 15-12-2016
% 
% version 12
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

pathToValidationFolder = '/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/04.VerificationAndValidation/04.DifferentVariables/nonRotatingMars';

currentVariableFolder = 'ErrorTolerance';

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
filePref = 0; % Choose a different file from the last one. So if this number is 1, you choose the penultimate file etc.


% Automatic latest file name (Source:
% https://www.youtube.com/watch?v=D8UkAOhsyDI)

% Write the paths to the folders
pathTSI = fullfile(pathToValidationFolder,'TSI',currentVariableFolder,currentCaseFolder);
pathRKF = fullfile(pathToValidationFolder,'RKF',currentVariableFolder,currentCaseFolder);

% Find all files ending with .csv in that folder
dTSI = dir(fullfile(pathTSI,'*.csv'));
dRKF = dir(fullfile(pathRKF,'*.csv'));

% Find the latest file
datesTSI = [dTSI.datenum]; % Get all the dates
[~,newestIndexTSI] = max(datesTSI); % Get the index of the newest date
newestFileTSI = dTSI(newestIndexTSI-filePref); % Get the data corresponding to the newest file
newestFileNameTSI = newestFileTSI.name; % Get the name of the newest file

datesRKF = [dRKF.datenum]; % Get all the dates
[~,newestIndexRKF] = max(datesRKF); % Get the index of the newest date
newestFileRKF = dRKF(newestIndexRKF-filePref); % Get the data corresponding to the newest file
newestFileNameRKF = newestFileRKF.name; % Get the name of the newest file



% Specify the paths
VariablePathTSI = fullfile(pathTSI,newestFileNameTSI); % Create the path for the TSI data

VariablePathRKF = fullfile(pathRKF,newestFileNameRKF); % Create the path for the RKF data


VariableVectorTSI = csvread(VariablePathTSI); % Read the file
VariableVectorRKF = csvread(VariablePathRKF); % Read the file


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
currentVariable = VariableVectorTSI(:,1); % Current variable values
Longitude = VariableVectorTSI(:,2); % Initial longitude [deg]
cpuTimeTSI = VariableVectorTSI(:,3); % The CPU time estimate for TSI [sec]
wallTimeTSI = VariableVectorTSI(:,4); % The wall time estimate for TSI [sec]
functionEvaluationsTSI = VariableVectorTSI(:,5); % The TSI function evaluations

% End states
xPositionTSI = VariableVectorTSI(:,7); % x-position TSI [km]
yPositionTSI = VariableVectorTSI(:,8); % y-position TSI [km]
zPositionTSI = VariableVectorTSI(:,9); % z-position TSI [km]
xVelocityTSI = VariableVectorTSI(:,10); % x-velocity TSI [km/s]
yVelocityTSI = VariableVectorTSI(:,11); % y-velocity TSI [km/s]
zVelocityTSI = VariableVectorTSI(:,12); % z-velocity TSI [km/s]
MassTSI = VariableVectorTSI(:,13); % Mass TSI [kg]


% Accuracy of the results and speed of convergence
TruthTSI = VariableVectorTSI(end,7:13);
% TruthTSI = VariableVectorRKF((find(VariableVectorRKF(:,1)==20)),7:13);

% Difference in metres w.r.t. the order 20 end state

xPositionTSIdifference = abs(TruthTSI(1)-xPositionTSI)*1000; % x-position TSI [m]
yPositionTSIdifference = abs(TruthTSI(2)-yPositionTSI)*1000; % y-position TSI [m]
zPositionTSIdifference = abs(TruthTSI(3)-zPositionTSI)*1000; % z-position TSI [m]
xVelocityTSIdifference = abs(TruthTSI(4)-xVelocityTSI)*1000; % x-velocity TSI [m/s]
yVelocityTSIdifference = abs(TruthTSI(5)-yVelocityTSI)*1000; % y-velocity TSI [m/s]
zVelocityTSIdifference = abs(TruthTSI(6)-zVelocityTSI)*1000; % z-velocity TSI [m/s]
MassTSIdifference = abs(TruthTSI(7)-MassTSI); % Mass TSI [kg]

xPositionTSIdifferenceRKF = abs(VariableVectorTSI(:,21))*1000; % x-position TSI [m]
yPositionTSIdifferenceRKF = abs(VariableVectorTSI(:,22))*1000; % y-position TSI [m]
zPositionTSIdifferenceRKF = abs(VariableVectorTSI(:,23))*1000; % z-position TSI [m]
xVelocityTSIdifferenceRKF = abs(VariableVectorTSI(:,24))*1000; % x-velocity TSI [m/s]
yVelocityTSIdifferenceRKF = abs(VariableVectorTSI(:,25))*1000; % y-velocity TSI [m/s]
zVelocityTSIdifferenceRKF = abs(VariableVectorTSI(:,26))*1000; % z-velocity TSI [m/s]
MassTSIdifferenceRKF = abs(VariableVectorTSI(:,27)); % Mass TSI [kg]

% Convergence through difference between each consecutive value

for k = 1:1:length(xPositionTSI)
   
    if k ~= length(xPositionTSI)
    xPositionTSIdifferenceCons(k) = abs(xPositionTSI(k+1)-xPositionTSI(k))*1000; % x-position TSI [m]
yPositionTSIdifferenceCons(k) = abs(yPositionTSI(k+1)-yPositionTSI(k))*1000; % y-position TSI [m]
zPositionTSIdifferenceCons(k) = abs(zPositionTSI(k+1)-zPositionTSI(k))*1000; % z-position TSI [m]
xVelocityTSIdifferenceCons(k) = abs(xVelocityTSI(k+1)-xVelocityTSI(k))*1000; % x-velocity TSI [m/s]
yVelocityTSIdifferenceCons(k) = abs(yVelocityTSI(k+1)-yVelocityTSI(k))*1000; % y-velocity TSI [m/s]
zVelocityTSIdifferenceCons(k) = abs(zVelocityTSI(k+1)-zVelocityTSI(k))*1000; % z-velocity TSI [m/s]
MassTSIdifferenceCons(k) = abs(MassTSI(k+1)-MassTSI(k)); % Mass TSI [kg]
    end

if k == length(xPositionTSI)
   
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
semilogx(currentVariable,cpuTimeTSI);

title(['CPU time TSI vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('CPU time TSI [sec]');

% axis([currentVariable(1) currentVariable(end) 0.05 0.1]) % First zoom
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom

% axis([currentVariable(1) currentVariable(end) 0.05 0.16]) % First zoom
% case 2
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom
% case 2


% legend('CPU time','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(2) % Wall time
semilogx(currentVariable,wallTimeTSI);

title(['Wall time TSI vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Wall time TSI [sec]');

% axis([currentVariable(1) currentVariable(end) 0.05 0.1]) % First zoom
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom

% axis([currentVariable(1) currentVariable(end) 0.05 0.16]) % First zoom
% case 2
% axis([currentVariable(1) currentVariable(end) 0.05 0.06]) % Second zoom
% case 2

% legend('Wall time','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(3) % Position
semilogx(currentVariable,xPositionTSI)

title(['x-Position vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('x-Position [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean(xPositionTSI); % Mean of the vector
STD = std(xPositionTSI);      % Standard deviation of the vector
dummyVector = xPositionTSI;    % Create a replacement vector that can be adjusted

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

axis([currentVariable(end) currentVariable(1) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom



% legend('x-Position','y-Position','z-Position','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(4)
semilogx(currentVariable,yPositionTSI)

title(['y-Position vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('y-Position [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean(yPositionTSI); % Mean of the vector
STD = std(yPositionTSI);      % Standard deviation of the vector
dummyVector = yPositionTSI;    % Create a replacement vector that can be adjusted

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

axis([currentVariable(end) currentVariable(1) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Position','y-Position','z-Position','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(5)
semilogx(currentVariable,zPositionTSI)

title(['z-Position vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('z-Position [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean(zPositionTSI); % Mean of the vector
STD = std(zPositionTSI);      % Standard deviation of the vector
dummyVector = zPositionTSI;    % Create a replacement vector that can be adjusted

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

axis([currentVariable(end) currentVariable(1) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Position','y-Position','z-Position','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(6) % Velocity
semilogx(currentVariable,xVelocityTSI)
% hold on
% plot(currentVariable,yVelocityTSI)
% hold on
% plot(currentVariable,zVelocityTSI)

title(['x-Velocity vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('x-Velocity [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean(xVelocityTSI); % Mean of the vector
STD = std(xVelocityTSI);      % Standard deviation of the vector
dummyVector = xVelocityTSI;    % Create a replacement vector that can be adjusted

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

axis([currentVariable(end) currentVariable(1) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Velocity','y-Velocity','z-Velocity','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(7)
semilogx(currentVariable,yVelocityTSI)

title(['y-Velocity vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('y-Velocity [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean(yVelocityTSI); % Mean of the vector
STD = std(yVelocityTSI);      % Standard deviation of the vector
dummyVector = yVelocityTSI;    % Create a replacement vector that can be adjusted

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

axis([currentVariable(end) currentVariable(1) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Velocity','y-Velocity','z-Velocity','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(8)
semilogx(currentVariable,zVelocityTSI)

title(['z-Velocity vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('z-Velocity [km]');

% Determine the mean and the standard deviation of the vector
Mean = mean(zVelocityTSI); % Mean of the vector
STD = std(zVelocityTSI);      % Standard deviation of the vector
dummyVector = zVelocityTSI;    % Create a replacement vector that can be adjusted

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

axis([currentVariable(end) currentVariable(1) Meanadj-3*STDadj Meanadj+3*STDadj]) % Zoom

% legend('x-Velocity','y-Velocity','z-Velocity','Location','NorthEastOutside'); % Add a legend in the top right corner


figure(9)
semilogx(currentVariable,functionEvaluationsTSI)

title(['Function evaluations vs ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Function evaluations [-]');

% Determine the mean and the standard deviation of the vector
Mean = mean(functionEvaluationsTSI); % Mean of the vector
STD = std(functionEvaluationsTSI);      % Standard deviation of the vector
dummyVector = functionEvaluationsTSI;    % Create a replacement vector that can be adjusted

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

% axis([currentVariable(end) currentVariable(1) 0 Meanadj+3*STDadj]) % Zoom
% axis([currentVariable(1) currentVariable(end) 0 100]) % Zoom



figure(10)
semilogx(functionEvaluationsTSI,cpuTimeTSI)

title(['Function evaluations vs CPU time TSI']); % Give the figure a title
xlabel('Function evaluations [-]');
ylabel('CPU time TSI [sec]'); % Label the different axes


% axis([min(functionEvaluationsTSI)-5 100 0.05 0.09]) % Zoom

%% Speed of convergence and accuracy

figure(11)
loglog(currentVariable,xPositionTSIdifference)
hold on
loglog(currentVariable,yPositionTSIdifference)
hold on
loglog(currentVariable,zPositionTSIdifference)
hold on
loglog(currentVariable,xVelocityTSIdifference)
hold on
loglog(currentVariable,yVelocityTSIdifference)
hold on
loglog(currentVariable,zVelocityTSIdifference)
hold on
loglog(currentVariable,MassTSIdifference)


title(['Absolute differences w.r.t. the nominal case for every ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Absolute difference [m, m/s, kg]');

legend('x-Position','y-Position','z-Position','x-Velocity','y-Velocity','z-Velocity','Mass','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(12)
loglog(currentVariable,xPositionTSIdifferenceRKF)
hold on
loglog(currentVariable,yPositionTSIdifferenceRKF)
hold on
loglog(currentVariable,zPositionTSIdifferenceRKF)
hold on
loglog(currentVariable,xVelocityTSIdifferenceRKF)
hold on
loglog(currentVariable,yVelocityTSIdifferenceRKF)
hold on
loglog(currentVariable,zVelocityTSIdifferenceRKF)
hold on
loglog(currentVariable,MassTSIdifferenceRKF)


title(['Absolute differences w.r.t. the RKF case for every ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Absolute difference [m, m/s, kg]');

legend('x-Position','y-Position','z-Position','x-Velocity','y-Velocity','z-Velocity','Mass','Location','NorthEastOutside'); % Add a legend in the top right corner



figure(13)
loglog(currentVariable,xPositionTSIdifferenceCons)
hold on
loglog(currentVariable,yPositionTSIdifferenceCons)
hold on
loglog(currentVariable,zPositionTSIdifferenceCons)
hold on
loglog(currentVariable,xVelocityTSIdifferenceCons)
hold on
loglog(currentVariable,yVelocityTSIdifferenceCons)
hold on
loglog(currentVariable,zVelocityTSIdifferenceCons)
hold on
loglog(currentVariable,MassTSIdifferenceCons)


title(['Consecutive differences between ',currentVariableFolder]); % Give the figure a title
xlabel([currentVariableFolder,' ',currentVariableUnit]); % Label the different axes
ylabel('Consecutive differences [m, m/s, kg]');

legend('x-Position','y-Position','z-Position','x-Velocity','y-Velocity','z-Velocity','Mass','Location','NorthEastOutside'); % Add a legend in the top right corner

toc




