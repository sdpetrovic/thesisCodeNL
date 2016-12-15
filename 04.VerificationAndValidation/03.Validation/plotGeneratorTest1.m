% This is script lets me plot the 3-D graphs of the generated trajectories.
% For this the Cartesian generated data is used (also in the spherical
% computation case)
%
% Stacha Petrovic 31-08-2016
% 
% version 1
% 
% Linux Redhat 7.2
%
%
% 
%
%

close all
clear all
clc

tic

%% Required initial data

marsRadius = 3396; % [km]
marsStandardGravitationalParameter = 4.2828314e4;  % mu_M [km^3/s^2]

% Desired orbit parameters
desiredOrbitSemiMajorAxis = marsRadius+390.0; % a [kmMOLA]
desiredOrbitEccentricity = 0.0; % e [-]
desiredOrbitInclination = deg2rad(45); % i [rad]
desiredRAAN = 1.30693423357019; % Omega [rad]
desiredArgumentOfPerigee = 0.0; % omega [rad] (3.92078961310322 could be used)
desiredTrueAnomaly = 0.0; % theta [rad]


%% Create the plot data for the desired orbit

stepSizeTheta = 0.1; % Resolution [deg]

% Generate empty vectors
desiredOrbitXposition = zeros(360/stepSizeTheta,1); % x-position vector
desiredOrbitYposition = zeros(360/stepSizeTheta,1); % y-position vector
desiredOrbitZposition = zeros(360/stepSizeTheta,1); % z-position vector

desiredOrbitXvelocity = zeros(360/stepSizeTheta,1); % x-position vector
desiredOrbitYvelocity = zeros(360/stepSizeTheta,1); % y-velocity vector
desiredOrbitZvelocity = zeros(360/stepSizeTheta,1); % z-velocity vector



j = 1; %Counter used for the vectors

% i = 0;

for i=0:stepSizeTheta:(360.0-stepSizeTheta)
    
    [R,V] = randv(desiredOrbitSemiMajorAxis,desiredOrbitEccentricity,desiredOrbitInclination,desiredRAAN,desiredArgumentOfPerigee,deg2rad(i),marsStandardGravitationalParameter); % Convert kepler elements to cartesian coordinates
    
    
    desiredOrbitXposition(j,1) = R(1); % store the x-position
    desiredOrbitYposition(j,1) = R(2); % store the y-position
    desiredOrbitZposition(j,1) = R(3); % store the z-position
    
    j = j+1;
    
end

%% Get the data

% Specify the paths

RKFpath = fullfile('01.Version2.0FirstValidationTest','backupRKFFileAtDateAndTime_2016-09-02_11:18:13.csv'); % Create the path for the RKF file
TSIpath = fullfile('01.Version2.0FirstValidationTest','backupSpherical(Cart)TSIFileAtDateAndTime_2016-09-02_11:18:12.csv'); % Create the path for the TSI file
% referenceDataPath = fullfile('02.Version1SecondValidationTest','trajectoryOutputDataCase7_3_2016_v33.txt'); % Create the path for the reference trajectory
% RKFThrustPath = fullfile('01.Version2.0FirstValidationTest','Thrust.csv'); % Create the path for the RKF thrust file


  % Read the files
% CartesianRKFdata = dlmread(RKFpath,',', 1, 1); % Read the file
% CartesianTSIdata_spher = dlmread(TSIpath,',', 1, 1); % Read the file

CartesianRKFdata = csvread(RKFpath); % Read the file
CartesianTSIdata_spher = csvread(TSIpath); % Read the file
% ReferenceData = dlmread(referenceDataPath); % Read the file
% ThrustData = csvread(RKFThrustPath); % Read the file


%% Generate the required vectors

% RKF
RKFxPosition = CartesianRKFdata(:,2); % RKF x-position
RKFyPosition = CartesianRKFdata(:,3); % RKF y-position
RKFzPosition = CartesianRKFdata(:,4); % RKF z-position
RKFxVelocity = CartesianRKFdata(:,5); % RKF x-velocity
RKFyVelocity = CartesianRKFdata(:,6); % RKF y-velocity
RKFzVelocity = CartesianRKFdata(:,7); % RKF z-velocity

% TSI
TSIxPosition = CartesianTSIdata_spher(:,2); % TSI x-position
TSIyPosition = CartesianTSIdata_spher(:,3); % TSI y-position
TSIzPosition = CartesianTSIdata_spher(:,4); % TSI z-position
TSIxVelocity = CartesianTSIdata_spher(:,5); % TSI x-velocity
TSIyVelocity = CartesianTSIdata_spher(:,6); % TSI y-velocity
TSIzVelocity = CartesianTSIdata_spher(:,7); % TSI z-velocity

% % Reference trajectory
% RefxPosition = ReferenceData(:,2)*1e-3; % Ref x-position
% RefyPosition = ReferenceData(:,3)*1e-3; % Ref y-position
% RefzPosition = ReferenceData(:,4)*1e-3; % Ref z-position
% RefxVelocity = ReferenceData(:,5)*1e-3; % Ref x-velocity
% RefyVelocity = ReferenceData(:,6)*1e-3; % Ref y-velocity
% RefzVelocity = ReferenceData(:,7)*1e-3; % Ref z-velocity

% % RKF Thrust
% k = 1;
% for l=1:1:length(RKFxPosition)
% 
% ThrustXacc(l,1) = -ThrustData(k,1)*1e12; % RKF x acc thrust
% ThrustYacc(l,1) = -ThrustData(k,2)*1e12; % RKF y acc thrust
% ThrustZacc(l,1) = -ThrustData(k,3)*1e12; % RKF z acc thrust
% 
% k = k+13;
% 
% end

%% Plots

[xMars, yMars, zMars] = sphere; % Generating a unit sphere

% Giving the unit sphere the radius of Mars [km]
xMars = xMars*marsRadius; 
yMars = yMars*marsRadius; 
zMars = zMars*marsRadius; 

figure(1)
Mars = surf(xMars,yMars,zMars); % Plot "Mars"
hold on
% set(Mars,'FaceColor',[0.8 0.2 0 ],'FaceAlpha',0.5); % Give Mars a nice orangy colour
set(Mars,'FaceColor',[1 1 1 ],'FaceAlpha',0.5); % Give Mars a nice white colour
hold on
plot3(RKFxPosition,RKFyPosition,RKFzPosition,'g'); % Plot the RKF trajectory
hold on
plot3(TSIxPosition,TSIyPosition,TSIzPosition,'b'); % Plot the TSI trajectory
hold on
plot3(desiredOrbitXposition,desiredOrbitYposition,desiredOrbitZposition,'k--'); % Plot the desired orbit
% hold on
% plot3(RefxPosition,RefyPosition,RefzPosition,'r--'); % Plot the reference trajectory 
% hold on
% quiver3(RKFxPosition,RKFyPosition,RKFzPosition,ThrustXacc,ThrustYacc,ThrustZacc,'r'); % Plot the RKF thrust
% view(180,0);
% view(74.5,0); % see the trajectory from the side
% view(74.5,90);  % see the trajectory from the top
% view(90,0); % as seen from the x-axis
% view(180,0); % as seen from the y-axis
view(0,90); % as seen from the z-axis
% view(0,-90); % as seen from the negative z-axis

axis([-3000 1000 500 3500 0 3000]);


title('3-D trajectory plot over Mars'); % Give the figure a title
xlabel('x-position [km]'); % Label the different axes
ylabel('y-position [km]');
zlabel('z-position [km]');
legend('Mars','RKF','TSI','Desired orbit','Location','NorthEastOutside'); % Add a legend in the top right corner


figure(2) % Just the trajectories
plot3(RKFxPosition,RKFyPosition,RKFzPosition,'g'); % Plot the RKF trajectory
hold on
plot3(TSIxPosition,TSIyPosition,TSIzPosition,'b'); % Plot the TSI trajectory
hold on
plot3(desiredOrbitXposition,desiredOrbitYposition,desiredOrbitZposition,'k--'); % Plot the desired orbit
% hold on
% plot3(RefxPosition,RefyPosition,RefzPosition,'r--'); % Plot the reference trajectory 

axis([-3000 1000 500 3500 0 3000]); % Set specific axes

title('3-D trajectory plot (trajectories only)'); % Give the figure a title
xlabel('x-position [km]'); % Label the different axes
ylabel('y-position [km]');
zlabel('z-position [km]');
legend('RKF','TSI','Desired orbit','Location','NorthEastOutside'); % Add a legend in the top right corner

% Create the same plots including the velocity 
figure(3)
Mars = surf(xMars,yMars,zMars); % Plot "Mars"
hold on
set(Mars,'FaceColor',[0.8 0.2 0 ],'FaceAlpha',0.5); % Give Mars a nice orangy colour
hold on
quiver3(RKFxPosition,RKFyPosition,RKFzPosition,RKFxVelocity,RKFyVelocity,RKFzVelocity,'g'); % Plot the RKF velocity
hold on
quiver3(TSIxPosition,TSIyPosition,TSIzPosition,TSIxVelocity,TSIyVelocity,TSIzVelocity,'b'); % Plot the TSI velocity
hold on
plot3(desiredOrbitXposition,desiredOrbitYposition,desiredOrbitZposition,'k--'); % Plot the desired orbit

title('3-D velocity plot over Mars'); % Give the figure a title
xlabel('x-position [km]'); % Label the different axes
ylabel('y-position [km]');
zlabel('z-position [km]');
legend('Mars','RKF','TSI','Desired orbit','Location','NorthEastOutside'); % Add a legend in the top right corner


figure(4) % Just the trajectories
% plot3(RKFxPosition,RKFyPosition,RKFzPosition,'y'); % Plot the RKF trajectory
% hold on
% plot3(TSIxPosition,TSIyPosition,TSIzPosition,'r'); % Plot the TSI trajectory
% hold on
quiver3(RKFxPosition,RKFyPosition,RKFzPosition,RKFxVelocity,RKFyVelocity,RKFzVelocity,'g'); % Plot the RKF velocity
hold on
quiver3(TSIxPosition,TSIyPosition,TSIzPosition,TSIxVelocity,TSIyVelocity,TSIzVelocity,'b'); % Plot the TSI velocity
hold on
plot3(desiredOrbitXposition,desiredOrbitYposition,desiredOrbitZposition,'k--'); % Plot the desired orbit

axis([-4000 4000 -4000 4000 -4000 4000]); % Set specific axes

title('3-D velocity plot (velocity only)'); % Give the figure a title
xlabel('x-position [km]'); % Label the different axes
ylabel('y-position [km]');
zlabel('z-position [km]');
legend('RKF','TSI','Desired orbit','Location','NorthEastOutside'); % Add a legend in the top right corner

figure(5)
plot3(desiredOrbitXposition,desiredOrbitYposition,desiredOrbitZposition,'k--'); % Plot the desired orbit

axis([-4000 4000 -4000 4000 -4000 4000]); % Set specific axes

title('3-D desired orbit plot (orbit only)'); % Give the figure a title
xlabel('x-position [km]'); % Label the different axes
ylabel('y-position [km]');
zlabel('z-position [km]');
legend('Desired orbit','Location','NorthEastOutside'); % Add a legend in the top right corner

% figure(6)
% quiver3(RKFxPosition,RKFyPosition,RKFzPosition,ThrustXacc,ThrustYacc,ThrustZacc,'r'); % Plot the RKF thrust
% 
% axis([-4000 4000 -4000 4000 -4000 4000]); % Set specific axes
% 
% title('3-D thrust plot (thrust only)'); % Give the figure a title
% xlabel('x-position [km]'); % Label the different axes
% ylabel('y-position [km]');
% zlabel('z-position [km]');
% legend('Thrust','Location','NorthEastOutside'); % Add a legend in the top right corner

toc

%%

% % FULL SCALE RUN     Launch site: Nili Fossae at -0.6 km w.r.t. MOLA geoid
% 
% launchLat = 21.0;       % Latitude of the chosen launch site
% launchLong = 74.5;      % Longitude of the chosen launch site
% 
% LatRange = 20;          % Range of latitudes to be computed around the launch site
% LongRange = LatRange;   % Range of longitudes to be computed around the launch site
% 
% initial_lat = launchLat-LatRange;        % Starting latitude in degrees
% initial_long = launchLong-LongRange;       % Starting longitude in degrees
% 
% stepSize = 5;          % Lat and long difference per step in degrees
% 
% final_lat = launchLat+LatRange;         % Final latitude in degrees
% final_long = launchLong+LongRange; % Final longitude in degrees
% 
% % Read txt into cell A
% 
% fid = fopen('test.txt','r');
% i = 1;
% tline = fgetl(fid);
% A{i} = tline;
% while ischar(tline)
%     i = i+1;
%     tline = fgetl(fid);
%     A{i} = tline;
% end
% fclose(fid);
% 
% % How many runs are still to come?
% count = ((final_lat-initial_lat+stepSize)/stepSize)*((final_long-initial_long+stepSize)/stepSize);
% 
% countInitial = count;
% 
% % disp(['Only ' num2str(count) ' runs still to go']);
% 
% 
% for lat = initial_lat:stepSize:final_lat
%    
%     % Set the latitude
%     
%     A{38} = sprintf('  FLAT     = %d',lat);
% 
%     
%     for long = initial_long:stepSize:final_long
%      
%         % Display how many runs are still to come
%         
%         disp(['Number of runs still to go: ' num2str(count)]);  
%         
%     % Set the longitude   
%     
%     A{39} = sprintf('  FLON     = %d',long);  
%     
%     % Write cell A into txt
%     fid = fopen('test.txt','w');
% 
%     for i = 1:numel(A)
%         if A{i+1} == -1
%             fprintf(fid, '%s', A{i});
%             break
%         else
%             fprintf(fid,'%s\n', A{i});
%         end
%     end
%     
% 
%     
%     
%     %% Run the program
%     
%         
%       fileExePath = '/home/stachap/Documents/Thesis/04._SNOPT_and_Mars-GRAM_2005_JPL/MarsGram2005-1.3/src/marsgram_M05 < nameFile.txt'; % Locate the program to be run and feed the input file at the same time
%     
%       % This final method seems to work. The key was to put the name of the
%       % file in the input .txt file that the system function needs!
%       
%     status = system(fileExePath);            % Run the program
% 
%   %% Read the files
% Density_height_density = dlmread('Density.txt',' ', 1, 1); % Read the file
% 
% display 'Density.txt has been read'; 
% 
% TPresHgt_temp_pres = dlmread('TPresHgt.txt',' ', 1, 1); % Read the file
% 
% display 'TPresHgt.txt has been read'; 
% 
% %% Get the proper data
% 
% % Height
% 
% i = 1;          % Row
% j = 1;          % Column
% 
% DIMENSION_vect = Density_height_density(:,1);
% DIMENSION = length(DIMENSION_vect);
%    
%     
%     LatitudeVectorData = zeros(DIMENSION,1);
%     LongitudeVectorData = zeros(DIMENSION,1);
%     
%     k = 1;
%     
%     while (k<=DIMENSION)
%         
%        LatitudeVectorData(k) = lat;
%        LongitudeVectorData(k) = long;
%        
%        k = k+1;
%         
%     end
%     
%     
% HeightVectorData = zeros(DIMENSION,1); 
% 
% while (i <= DIMENSION)
%    
%     if Density_height_density(i,j)==0
%         
%         j = j+1;
%         
%     else
%         
%         HeightVectorData(i) = Density_height_density(i,j);
%         
%         i = i+1;
%     end
%     
%     if j == 4
%         
%         disp('whoops');
%         
%         break
%         
%     end
%     
% end
% 
% % display 'Height data has successfully been aquired'; 
% 
% % Density
% 
% i = 1;          % Row
% j = 1;          % Column
% 
% DIMENSION_vect = Density_height_density(:,1);
% DIMENSION = length(DIMENSION_vect);
% 
% DensityVectorData = zeros(DIMENSION,1); 
% 
% numCount = 0;
% numLim = 2;
% 
% while (i <= DIMENSION)
%     
%        if Density_height_density(i,j)==0
%         
%         j = j+1;
%         
%        elseif Density_height_density(i,j)~=0 && numCount<numLim
%            
%            numCount = numCount+1;
%            
%            j = j+1;
%            
%        else
%         
%         DensityVectorData(i) = Density_height_density(i,j);
%         
%         numCount = 0;
%         
%         j = 1;
%         
%         i = i+1;
%         
%     end
%     
% %     if j == 12
% %         
% %         disp('whoops');
% %         
% %         break
% %         
% %     end 
%     
%     
%     
% end
% 
% 
% 
% % display 'Density data has successfully been aquired'; 
% 
% % Temperature
% 
% i = 1;          % Row
% j = 1;          % Column
% 
% DIMENSION_vect = TPresHgt_temp_pres(:,1);
% DIMENSION = length(DIMENSION_vect);
% 
% TemperatureVectorData = zeros(DIMENSION,1); 
% numCount = 0;
% numLim = 1;
% 
% while (i <= DIMENSION)
%     
%        if TPresHgt_temp_pres(i,j)==0
%         
%         j = j+1;
%         
%        elseif TPresHgt_temp_pres(i,j)~=0 && numCount<numLim
%            
%            numCount = numCount+1;
%            
%            j = j+1;
%            
%        else
%         
%         TemperatureVectorData(i) = TPresHgt_temp_pres(i,j);
%         
%         numCount = 0;
%         
%         j = 1;
%         
%         i = i+1;
%         
%     end
%     
% %     if j == 11
% %         
% %         disp('whoops');
% %         
% %         break
% %         
% %     end 
%     
%     
% end
% 
% 
% 
% % display 'Temperature data has successfully been aquired'; 
% 
% % Pressure
% 
% i = 1;          % Row
% j = 1;          % Column
% 
% DIMENSION_vect = TPresHgt_temp_pres(:,1);
% DIMENSION = length(DIMENSION_vect);
% 
% PressureVectorData = zeros(DIMENSION,1); 
% 
% numCount = 0;
% numLim = 2;
% 
% while (i <= DIMENSION)
%     
%        if TPresHgt_temp_pres(i,j)==0
%         
%         j = j+1;
%         
%        elseif TPresHgt_temp_pres(i,j)~=0 && numCount<numLim
%            
%            
%            numCount = numCount+1;
%            
%            j = j+1;
%            
%        else
%         
%         PressureVectorData(i) = TPresHgt_temp_pres(i,j);
%         
%         numCount = 0;
%         
%         j = 1;
%         
%         i = i+1;
%         
%     end
%     
% %     if j == 13
% %         
% %         disp('whoops');
% %         
% %         break
% %         
% %     end 
%     
%     
%     
% end
% 
% 
% % display 'Pressure data has successfully been aquired'; 
% 
% %% Put the data table together
% 
% % Altitude, Latitude, Longitude, Temperature, Pressure and Density
% % respectively
% 
% H_T_P_rho_dataMatrix = zeros(DIMENSION,6);
% 
% H_T_P_rho_dataMatrix(:,1) = HeightVectorData;
% H_T_P_rho_dataMatrix(:,2) = LatitudeVectorData;
% H_T_P_rho_dataMatrix(:,3) = LongitudeVectorData;
% H_T_P_rho_dataMatrix(:,4) = TemperatureVectorData;
% H_T_P_rho_dataMatrix(:,5) = PressureVectorData;
% H_T_P_rho_dataMatrix(:,6) = DensityVectorData;
% 
% % display 'Data matrix has been set-up';
% 
% %% Write .csv file
% 
% 
% 
% pathname = fullfile('TestTableHTPrho_fullScale_8_fiveDegreesAccurate_launchSitePM20deg.csv');
% 
% 
% dlmwrite(pathname,H_T_P_rho_dataMatrix, '-append');   % Used to append to
% %the first file
% 
% display 'Data has been added to file';
%    
% 
% 
%     count = count -1;
%     
%     
%         
%     end
%     
% end
% 
% display 'Done!! You now have a table :D'


