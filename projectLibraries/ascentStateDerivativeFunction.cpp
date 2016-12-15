/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      160517    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

#include "ascentStateDerivativeFunction.h"

/// The state derivative function ///
/// \brief ascentStateDerivativeFunction The function that determines the state derivatives
/// \param planet_          The celestial body class, in this case Mars
/// \param MAV_             The Mars ascent vehicle class
/// \param currentStateAndTime_ The state and time class
/// \return
///

const tudat::basic_mathematics::Vector7d ascentStateDerivativeFunction(const celestialBody& planet_, const MarsAscentVehicle& MAV_, const StateAndTime& currentStateAndTime_){



    // Create variables to be used in this function

        celestialBody Mars = planet_;                               // The celestial body class
        MarsAscentVehicle MAV = MAV_;                                 // The MAV class
        StateAndTime stateAndTime = currentStateAndTime_;      // The current state and time class

        const double rotationalVelocityMars = Mars.rotationalVelocity();
        const double primeMeridianAngle = Mars.primeMeridianAngle();
        const double inertialFrameTime = Mars.inertialFrameTime();

/// Compute current spherical state ///

    tudat::basic_mathematics::Vector7d currentStateAndTime = stateAndTime.getCurrentState();

    const double xPosition = currentStateAndTime(0);            // x position coordinate definition
    const double yPosition = currentStateAndTime(1);            // y position coordinate definition
    const double zPosition = currentStateAndTime(2);            // z position coordinate definition
    const double xVelocity = currentStateAndTime(3);            // x velocity coordinate definition
    const double yVelocity = currentStateAndTime(4);            // y velocity coordinate definition
    const double zVelocity = currentStateAndTime(5);            // z velocity coordinate definition
    const double massMAV = currentStateAndTime(6);              // MAV mass definition
    const double currentTime = stateAndTime.getCurrentTime();   // current time definition

/*    /// Setting the data collection file for RKF and inserting the first values ///

        // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
        const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/04.Verification/00.RKFintermediateStateOutput/";


        // Set output format for matrix output.
        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

        // Set absolute path to file containing the data.
        std::string dataAbsolutePath = outputDirectory + "test1FullIntegrationRK4stateAndTime.csv";

        // Create a row vector for the storing of the data
        Eigen::MatrixXd outputVector = Eigen::MatrixXd::Zero(1,8); // Create a row vector for the storing of the data

        // Filling the output vector
        outputVector(0,0) = currentTime;   // Storing the initial time
        outputVector(0,1) = xPosition;   // Storing the initial x position
        outputVector(0,2) = yPosition;   // Storing the initial y position
        outputVector(0,3) = zPosition;   // Storing the initial z position
        outputVector(0,4) = xVelocity;   // Storing the initial x velocity
        outputVector(0,5) = yVelocity;   // Storing the initial y velocity
        outputVector(0,6) = zVelocity;   // Storing the initial z velocity
        outputVector(0,7) = massMAV;   // Storing the initial MAV mass


        // Storing the data

        std::ifstream ifile(dataAbsolutePath.c_str()); // Check it as an input file

        bool fexists = false;   // Set the default to "It does not exist"

        if (ifile){         // Attempt to open the file


           fexists = true;      // If the file can be opened it must exist

           ifile.close();   // Close the file

        }


        // If so: error and create temporary file, if not: create new file and put data in

        if (fexists == true){


            // Get time //

            time_t rawtime;
            struct tm * timeinfo;

            time ( &rawtime );
            timeinfo = localtime ( &rawtime );
    //        printf ( "Current local time and date: %s", asctime (timeinfo) );         // (from stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c and from stackoverflow.com/questions/1442116/how-to-get-date-and-time-value-in-c-program)
    //        std::cout<<"time = "<<time ( &rawtime )<<std::endl;
    //        std::cout<< (timeinfo->tm_year+1900) << "-"
    //                 << (timeinfo->tm_mon + 1) << "-"
    //                 << timeinfo->tm_mday << " "
    //                 << timeinfo->tm_hour << ":"
    //                 << timeinfo->tm_min << ":"
    //                 << timeinfo->tm_sec <<std::endl;

            ostringstream ConvertYear;
            ConvertYear << (timeinfo->tm_year+1900);
            std::string currentYear = ConvertYear.str();

            ostringstream ConvertMonth;
            ConvertMonth << (timeinfo->tm_mon+1);
            std::string currentMonth = ConvertMonth.str();

            ostringstream ConvertDay;
            ConvertDay << timeinfo->tm_mday;
            std::string currentDay = ConvertDay.str();

            ostringstream ConvertHour;
            ConvertHour << timeinfo->tm_hour;
            std::string currentHour = ConvertHour.str();

    //        std::cout<<"currentHour = "<<currentHour<<std::endl;

            ostringstream ConvertMin;
            ConvertMin << timeinfo->tm_min;
            std::string currentMin = ConvertMin.str();

    //        std::cout<<"currentMin = "<<currentMin<<std::endl;

            ostringstream ConvertSec;
            ConvertSec << timeinfo->tm_sec;
            std::string currentSec = ConvertSec.str();


            // Making sure each one of the representation has at two numbers
            if(currentMonth.size() == 1){
                currentMonth = "0" + currentMonth;
            }


            if(currentDay.size() == 1){
                currentDay = "0" + currentDay;
            }

            if(currentSec.size() == 1){
                currentSec = "0" + currentSec;
            }

            if (currentMin.size() == 1){
                currentMin = "0" + currentMin;
            }

            if (currentHour.size() == 1){
                currentHour = "0" + currentHour;
            }


    //        std::cout<<"The length of currentSec = "<<currentSec.size()<<" and the value = "<<currentSec<<std::endl;

            // Create new file name

            std::string ComputerTimeString = currentYear + "-" + currentMonth + "-" + currentDay + "_" + currentHour + ":" + currentMin + ":" + currentSec;  // Convert to string and store

            std::string newFileName = "backupRKFFileAtDateAndTime_" + ComputerTimeString + ".csv";

        std::cerr<<"The file name that you have chosen already exists, a new file with name "<<newFileName<<" will be created to store the data for now"<<std::endl;

        // Set new absolute path to file containing the data.
        dataAbsolutePath = outputDirectory + newFileName;

        // Export the data.
        std::ofstream exportFile1( dataAbsolutePath.c_str( ) ); // Make the new file
        std::cout<<"New file called "<<dataAbsolutePath<<" has been created"<<std::endl;
        exportFile1 << outputVector.format( csvFormat );          // Store the new values
        exportFile1.close( );   // Close the file


    }
            else{

            // Export the data.
            std::ofstream exportFile1( dataAbsolutePath.c_str( ) ); // Make the new file
            std::cout<<"New file called "<<dataAbsolutePath<<" has been created"<<std::endl;
            exportFile1 << outputVector.format( csvFormat );          // Store the new values
            exportFile1.close( );   // Close the file
        };
        //*/


// Computations

    const double Radius = sqrt(xPosition*xPosition+yPosition*yPosition+zPosition*zPosition);         // r [km]

    const double inertialVelocity = sqrt(xVelocity*xVelocity+yVelocity*yVelocity+zVelocity*zVelocity);       // V_I [km/s]

    const double inertialLongitude = atan2(yPosition,xPosition);         // lambda [rad]

    const double Latitude = asin(zPosition/Radius);              // delta [rad]

    const double rotationalLongitude = inertialLongitude-rotationalVelocityMars*(inertialFrameTime+currentTime)+primeMeridianAngle;  // tau [rad]

    double inertialLongitudeChange_;       // lambda_dot [rad/s] (placeholder)

    // Avoiding singularities
    if ((xPosition*xPosition+yPosition*yPosition) == 0){

        inertialLongitudeChange_ = 0;
    }
    else {
        inertialLongitudeChange_ = (xPosition*yVelocity-yPosition*xVelocity)/(xPosition*xPosition+yPosition*yPosition);
    };

    const double inertialLongitudeChange = inertialLongitudeChange_;        // lambda_dot [rad/s] (actual parameter)

    double rotationalLongitudeChange_ = inertialLongitudeChange-rotationalVelocityMars;     // tau_dot [rad/s] (placeholder)

    if (rotationalLongitudeChange_<=1e-15){  // Setting the accuracy to 1e-15 to avoid problems in the beginning with rounding errors...

        rotationalLongitudeChange_ = 0;

    };

    const double rotationalLongitudeChange = rotationalLongitudeChange_;    // tau_dot [rad/s] (actual parameter)

/*    /// Debug ///

    std::cout<<"rotationalVelocityMars = "<<rotationalVelocityMars<<std::endl;
    std::cout<<"inertialLongitudeChange = "<<inertialLongitudeChange<<std::endl;
    std::cout<<"rotationalVelocityMars-7.088e-05 = "<<rotationalVelocityMars-7.088e-05<<std::endl;
    std::cout<<"inertialLongitudeChange-7.088e-05 = "<<inertialLongitudeChange-7.088e-05<<std::endl;
    std::cout<<"(xPosition*yVelocity-yPosition*xVelocity) = "<<(xPosition*yVelocity-yPosition*xVelocity)<<std::endl;
    std::cout<<"(xPosition*xPosition+yPosition*yPosition) = "<<(xPosition*xPosition+yPosition*yPosition)<<std::endl;
//*/


    const double RadiusChange = (xPosition*xVelocity+yPosition*yVelocity+zPosition*zVelocity)/(Radius);      // radial velocity [km/s]

    double LatitudeChange_; // delta_dot [rad/s] (placeholder)

    if ((Radius*Radius*sqrt(1-(zPosition/Radius)*(zPosition/Radius))) == 0){
        LatitudeChange_ = 0;
    }
    else{
        LatitudeChange_ = (Radius*zVelocity-zPosition*RadiusChange)/(Radius*Radius*sqrt(1-(zPosition/Radius)*(zPosition/Radius)));
    };

    const double LatitudeChange = LatitudeChange_;  // delta_dot [rad/s] (actual parameter)

    // Avoid cosine round-off errors
    double localMarsRotationalVelocity = rotationalVelocityMars*Radius*cos(Latitude);  // V_M [km/s]

    if (abs(cos(Latitude))<6.2e-17){
      localMarsRotationalVelocity = 0;
    }

//    std::cout<<"RadiusChange = "<<RadiusChange<<std::endl;
//    std::cout<<"inertialVelocity = "<<inertialVelocity<<std::endl;

    double inertialFlightPathAngle = asin(RadiusChange/inertialVelocity);          // gamma_I [rad]
    // Avoid singularities
    if (inertialVelocity == 0){
        inertialFlightPathAngle = 0;
    }

    const double inertialAzimuth = atan2((inertialLongitudeChange*cos(Latitude)),LatitudeChange);    // chi_I [rad]

    const double rotationalVelocity = sqrt(localMarsRotationalVelocity*localMarsRotationalVelocity+inertialVelocity*inertialVelocity-2*localMarsRotationalVelocity*inertialVelocity*cos(inertialFlightPathAngle)*sin(inertialAzimuth));  // V_R [km/s]

    double rotationalFlightPathAngle_; // gamma_R [rad]  (placeholder)

    if (rotationalVelocity == 0){       // Setting the initial flight path angle in the rotational frame to 90 deg (or pi/2)

        rotationalFlightPathAngle_ = tudat::mathematical_constants::LONG_PI/2;
    }
    else {
        rotationalFlightPathAngle_ = asin(RadiusChange/rotationalVelocity);
    };

    const double rotationalFlightPathAngle = rotationalFlightPathAngle_;    // gamma_R [rad] (actual parameter)

    /// Debug ///
//    std::cout<<"rotationalFlightPathAngle = "<<rotationalFlightPathAngle<<std::endl;
//    std::cout<<"rotationalVelocity = "<<rotationalVelocity<<std::endl;
    /// Debug ///

    const double rotationalAzimuth = atan2((rotationalLongitudeChange*cos(Latitude)),LatitudeChange);    // chi_R [rad]

    // Check output
    std::cerr<<"Radius = "<<Radius<<std::endl;
    std::cout<<"inertialVelocity = "<<inertialVelocity<<std::endl;
    std::cout<<"inertialLongitude = "<<inertialLongitude<<std::endl;
    std::cout<<"Latitude = "<<Latitude<<std::endl;
    std::cout<<"rotationalLongitude = "<<rotationalLongitude<<std::endl;
    std::cout<<"inertialLongitudeChange = "<<inertialLongitudeChange<<std::endl;
    std::cout<<"rotationalLongitudeChange = "<<rotationalLongitudeChange<<std::endl;
    std::cout<<"RadiusChange = "<<RadiusChange<<std::endl;
    std::cout<<"LatitudeChange = "<<LatitudeChange<<std::endl;
    std::cout<<"localMarsRotationalVelocity = "<<localMarsRotationalVelocity<<std::endl;
    std::cout<<"inertialFlightPathAngle = "<<inertialFlightPathAngle<<std::endl;
    std::cout<<"inertialAzimuth = "<<inertialAzimuth<<std::endl;
    std::cout<<"rotationalVelocity = "<<rotationalVelocity<<std::endl;
    std::cout<<"rotationalFlightPathAngle = "<<rotationalFlightPathAngle<<std::endl;
    std::cout<<"rotationalAzimuth = "<<rotationalAzimuth<<std::endl;

//*/


/// Testing the local air temperature function ///

    const double currentAltitude = Radius-Mars.bodyReferenceRadius();

    const double currentTemperature = air_temperature::airTemperature(Mars.temperaturePolyCoefficients(), Mars.temperatureAltitudeRanges(),currentAltitude);

/*   // Check output
    std::cout<<"currentAltitude = "<<currentAltitude<<std::endl;
    std::cout<<"currentTemperature = "<<currentTemperature<<std::endl;
    std::cout<<"Radius = "<<Radius<<std::endl;
    std::cout<<"R_MOLA = "<<Mars.bodyReferenceRadius()<<std::endl;
    std::cout<<"Radius-3395.4 = "<<Radius-3395.4<<std::endl;
    std::cout<<"R_MOLA-3396 = "<<Mars.bodyReferenceRadius()-3396<<std::endl;
//*/

/// Testing the local air density function ///

    const double currentDensity= air_density::airDensity(Mars.densityPolyCoefficients(),  currentAltitude);

//    std::cout<<"The current air density = "<<currentDensity<<std::endl;

/// Testing the ascentDragForce function ///

/*    // Initially only with machNumber
    const double machNumber = Drag::ascentDragForce(rotationalVelocity,currentTemperature,Mars.adiabeticIndex(),Mars.specificGasConstant());

    std::cout<<"The current Mach number = "<<machNumber<<std::endl;


/// Testing the dragCoefficient function ///

    const double currentDragCoefficient = Drag::dragCoefficient(machNumber,MAV.dragCoefficientPolyCoefficients(),MAV.dragCoefficientMachRanges());

    std::cout<<"The current drag coefficient = "<<currentDragCoefficient<<std::endl;

    // Drag coefficient function as part of the drag function
    const double currentDragCoefficient = Drag::ascentDragForce(rotationalVelocity,currentTemperature,Mars.adiabeticIndex(),Mars.specificGasConstant(),
                                                                MAV.dragCoefficientPolyCoefficients(),MAV.dragCoefficientMachRanges());


    std::cout<<"The current drag coefficient = "<<currentDragCoefficient<<std::endl;
//*/

    const double currentDrag = Drag::ascentDragForce(rotationalVelocity,currentTemperature,Mars.adiabeticIndex(),Mars.specificGasConstant(),
                                                     MAV.dragCoefficientPolyCoefficients(),MAV.dragCoefficientMachRanges(),
                                                     MAV.referenceArea(), currentDensity);

//    std::cout<<"currentDrag = "<<currentDrag<<std::endl;



    /// Thrust acceleration in B-frame ///   thrustAccelerationsBframe

    const Eigen::Vector3d thrustAccelerationsPframe = Eigen::Vector3d((MAV.Thrust()/massMAV),0,0);            // THIS HAS TO BE CHANGED IN THE FUTURE TO INCLUDE A WIDE RANGE OF THRUST AZIMUTH AND ELEVATION ANGLES!!!

    const double thrustAzimuthTestDeg = 0;             // thrust azimuth gimbal angle [Deg] 10 for testing
    const double thrustElevationTestDeg = 0;            // thrust elevation gimbal angle [Deg] 5 for testing

    const double thrustAzimuthTest = deg2rad(thrustAzimuthTestDeg);     // thrust azimuth gimbal angle [rad]
    const double thrustElevationTest = deg2rad(thrustElevationTestDeg); // thrust elevation gimbal angle [rad]


    const Eigen::Vector3d thrustAccelerationsBframe = getPropulsionToBodyFrameTransformationMatrix(thrustAzimuthTest,thrustElevationTest)*thrustAccelerationsPframe;

//    std::cout<<"The thrust accelerations in the B-frame are "<<thrustAccelerationsBframe<<std::endl;

    /// Drag acceleration in B-frame ///

    const Eigen::Vector3d dragAccelerationsBframe = Eigen::Vector3d((-currentDrag/massMAV),0,0);

//    std::cout<<"The drag accelerations in the B-frame are "<<dragAccelerationsBframe<<std::endl;

/// Transfer to the inertial frame ///
///
    /// Accelerations in the V-frame ///

    // Thrust accelerations from B-frame to V-frame
    const Eigen::Vector3d thrustAccelerationsVframe = tudat::reference_frames::getTrajectoryToLocalVerticalFrameTransformationMatrix(rotationalFlightPathAngle,rotationalAzimuth)*thrustAccelerationsBframe;

    // Drag accelerations from B-frame to V-frame
    const Eigen::Vector3d dragAccelerationsVframe = tudat::reference_frames::getTrajectoryToLocalVerticalFrameTransformationMatrix(rotationalFlightPathAngle,rotationalAzimuth)*dragAccelerationsBframe;

//    std::cout<<"The thrust accelerations in the V-frame are "<<thrustAccelerationsVframe<<std::endl;
//    std::cout<<"The drag accelerations in the V-frame are "<<dragAccelerationsVframe<<std::endl;


    /// Accelerations in the R-frame ///

    // Thrust accelerations from V-frame to R-frame
    const Eigen::Vector3d thrustAccelerationsRframe = tudat::reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(rotationalLongitude,Latitude)*thrustAccelerationsVframe;

    // Drag accelerations from V-frame to R-frame
    const Eigen::Vector3d dragAccelerationsRframe = tudat::reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(rotationalLongitude,Latitude)*dragAccelerationsVframe;

//    std::cout<<"The thrust accelerations in the R-frame are "<<thrustAccelerationsRframe<<std::endl;
//    std::cout<<"The drag accelerations in the R-frame are "<<dragAccelerationsRframe<<std::endl;

    /// Accelerations in the I-frame ///

    const double angleItoR = rotationalVelocityMars*(inertialFrameTime+currentTime)-primeMeridianAngle;

    // Thrust accelerations from R-frame to I-frame
    const Eigen::Vector3d thrustAccelerationsIframe = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(angleItoR)*thrustAccelerationsRframe;

    // Drag accelerations from R-frame to I-frame
    const Eigen::Vector3d dragAccelerationsIframe = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(angleItoR)*dragAccelerationsRframe;

//    std::cout<<"The thrust accelerations in the I-frame are "<<thrustAccelerationsIframe<<std::endl;
//    std::cout<<"The drag accelerations in the I-frame are "<<dragAccelerationsIframe<<std::endl;

//    std::cout<<"The thrust+drag accelerations in the I-frame are "<<thrustAccelerationsIframe+dragAccelerationsIframe<<std::endl;


/// Compute gravitational acceleration ///

    Eigen::Vector3d gravAccelerationsIframe_;        // Define the placeholder gravitational acceleration vector

    for (int i = 0; i<3;i++){

        gravAccelerationsIframe_(i)=-Mars.standardGravitationalParameter()*(currentStateAndTime(i)/(pow(Radius,3)));

//        std::cout<<"gravAccelerationsIframe_("<<i<<") = "<<gravAccelerationsIframe_(i)<<std::endl;
    };

    const Eigen::Vector3d gravAccelerationsIframe = gravAccelerationsIframe_;         // Actual gravitational acceleration vector


//    std::cout<<"The gravitational accelerations in the I-frame are "<<gravAccelerationsIframe<<std::endl;

/// Compute total acceleration ///

    const Eigen::Vector3d totalAccelerationsIframe = gravAccelerationsIframe+dragAccelerationsIframe+thrustAccelerationsIframe;



/// Compute the mass flow rate ///

    const double massFlowRate = -MAV.Thrust()/(tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION*MAV.specificImpulse());



//    std::cout<<"The velocities in the I-frame are "<<stateAndTime.getCurrentVelocity()<<std::endl;
//    std::cout<<"The total accelerations in the I-frame are "<<totalAccelerationsIframe<<std::endl;
//    std::cout<<"The mass flow rate = "<<massFlowRate<<std::endl;


///// Define the output vector and fill it ///


    tudat::basic_mathematics::Vector7d stateDerivativeVector;

    stateDerivativeVector(0) = xVelocity;
    stateDerivativeVector(1) = yVelocity;
    stateDerivativeVector(2) = zVelocity;
    stateDerivativeVector(3) = totalAccelerationsIframe(0);
    stateDerivativeVector(4) = totalAccelerationsIframe(1);
    stateDerivativeVector(5) = totalAccelerationsIframe(2);
    stateDerivativeVector(6) = massFlowRate;


/*    /// Storing the values to the file ///


//    // Check if the file already exists.


//    std::ifstream ifile2(dataAbsolutePath.c_str()); // Check it as an input file

//    fexists = false;   // Set the default to "It does not exist"

//    if (ifile2){         // Attempt to open the file


//       fexists = true;      // If the file can be opened it must exist

//       ifile2.close();   // Close the file

//    }


//    // If so: append, if not: create new file and put data in

//    if (fexists == true){

//        // Export the Taylor Series Coefficients matrix.
//        std::ofstream exportFile1;                          // Define the file as an output file


//        exportFile1.open(dataAbsolutePath.c_str(),std::ios_base::app);      // Open the file in append mode

//        exportFile1 << "\n";                                            // Make sure the new matrix start on a new line

//        exportFile1 << stateDerivativeVector.format( csvFormat ); // Add the new values

//        std::cout<<"The file called "<<dataAbsolutePath<<" has been appended"<<std::endl;


//        exportFile1.close( );   // Close the file
//}
//        else{

//        std::cerr<<"Error: values could not be stored because storage file does not exist"<<std::endl;
//    };
    //*/



    return stateDerivativeVector;
    } // end of the state derivative function
