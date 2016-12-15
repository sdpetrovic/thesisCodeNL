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
 *      160427    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

// This is a test main file to test the different class files, header/source files to see if any output is produced (verification)

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <stdio.h>

#include <time.h>   // To determine the current computer time
#include <sys/time.h> // To determine the current computer time

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <cmath>

//#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>
#include <Tudat/Mathematics/BasicMathematics/coordinateConversions.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h>
//#include <tudatApplications/thesisProject/referenceFrameTransformationsUpdated.h>

/// Testing the celestial body class ///

//#include <thesisProject/celestialbody.h>          // Original test file
//#include <thesisProject/celestialBodyTestOne.h>     // First test
//#include <thesisProject/celestialBodyTestTwo.h>     // Second test
//#include <thesisProject/celestialBodyTestThree.h>   // Third test
//#include <thesisProject/celestialBodyTestFour.h>   // Fourth test

#include <thesisProject/celestialBody.h>            // Final version

/// Testing the vehicle class ///

//#include <thesisProject/MarsAscentVehicleTestOne.h> // First Test
//#include <thesisProject/MarsAscentVehicleTestTwo.h> // Second Test

#include <thesisProject/MarsAscentVehicle.h>    // Final version

/// Testing the current state and time and its updater ///

//#include <thesisProject/StateAndTime.h>           // Original test file
//#include <thesisProject/stateAndTimeTestTwo.h>      // Second test file
//#include <thesisProject/stateAndTimeTestThree.h>    // Third test file

#include <thesisProject/stateAndTime.h>             // Final version

/// Testing the auxiliary equations ///
#include <thesisProject/Auxiliary.h>                // Original test file


/// Testing the basic recurrence relations ///
#include <thesisProject/projectLibraries/basicRecurrenceRelations.h>               // Original test file

/// Testing all recurrence relations ///
#include <thesisProject/projectLibraries/allRecurrenceRelations.h>          // Original test file

/// Testing the stepSize class ///
#include <thesisProject/StepSize.h>             // Original test file

/// Testing the actual Taylor Series integration fucntion ///
#include <thesisProject/projectLibraries/TaylorSeriesIntegration.h>             // Original test file

/// Testing the other required functions ///
#include <thesisProject/projectLibraries/otherRequiredFunctions.h>              // Original test file


// testing




int main()

{

std::cout<<setprecision(15)<<"Setting output precision to 15"<<std::endl;

    /// Setting the Celestial Body class ///


//    // First test

//    celestialBody Mars;

    // Second test

//    const std::string planet = "Mars";
//    const std::string planet = "Venus";

    celestialBody Mars;

//    const double adiabeticIndex = Mars.adiabeticIndex();
//    const double specificGasConstant = Mars.specificGasConstant();
//    const double standardGravitationalParameter = Mars.standardGravitationalParameter();
    const double rotationalVelocity = Mars.rotationalVelocity();
    const double primeMeridianAngle = Mars.primeMeridianAngle();
    const double inertialFrameTime = Mars.inertialFrameTime();

   const double bodyReferenceRadius = Mars.bodyReferenceRadius();

//    celestialBody Mars(planet);
//    Mars.setPlanet(planet);  // Does not exist in the class anymore!


    /// Setting the vehicle class ///

    MarsAscentVehicle MAV;



  /// Initial conditions ///

    // Launch site characteristics

//    const double initialAltitude = -0.6e3;             // Starting altitude [m MOLA]
    const double initialAltitude = -0.6;                 // Starting altitude [km MOLA]
    const double initialLatitudeDeg = 21;               // Starting latitude [deg]
    const double initialLongitudeDeg = 74.5;            // Starting longitude [deg]

//    const double initialLatitude = initialLatitudeDeg*tudat::mathematical_constants::LONG_PI/180;       // Starting latitude [rad]
//    const double initialLongitude = initialLongitudeDeg*tudat::mathematical_constants::LONG_PI/180;     // Starting longitude [rad]

    const double initialLatitude = deg2rad(initialLatitudeDeg);       // Starting latitude [rad]
    const double initialLongitude = deg2rad(initialLongitudeDeg);     // Starting longitude [rad]

    const double initialRadius = bodyReferenceRadius+initialAltitude;               // Starting radius in km
//        const double initialRadius = bodyReferenceRadius+initialAltitude;               // Starting radius in m

        // Converting the initial spherical position to cartesian position using the standard convertSphericalToCartesian function of Tudat
        // Please note that this function requires the zenith angle as input which is pi/2-latitude!

    Eigen::Vector3d initialCartesianPositionRotationalFrame = Eigen::Vector3d::Zero(3);

      initialCartesianPositionRotationalFrame(0) = initialRadius*cos(initialLatitude)*cos(initialLongitude); // x_R
      initialCartesianPositionRotationalFrame(1) = initialRadius*cos(initialLatitude)*sin(initialLongitude); // y_R
      initialCartesianPositionRotationalFrame(2) = initialRadius*sin(initialLatitude); // z_R

    const Eigen::Vector3d initialCartesianPositionInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocity*inertialFrameTime-primeMeridianAngle)*initialCartesianPositionRotationalFrame;


    // Compute initial velocity in y-direction as seen from the launch site in the inertial frame

    const Eigen::Vector3d initialVelocityLaunchSite = Eigen::Vector3d(0,(rotationalVelocity*initialRadius*cos(initialLatitude)),0);

    const Eigen::Vector3d initialVelocityInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocity*inertialFrameTime-primeMeridianAngle+initialLongitude)*initialVelocityLaunchSite;

    /// Setting StateAndTime class using modified vector ///

    tudat::basic_mathematics::Vector7d aState;

    aState(0) = initialCartesianPositionInertialFrame(0);
    aState(1) = initialCartesianPositionInertialFrame(1);
    aState(2) = initialCartesianPositionInertialFrame(2);
    aState(3) = initialVelocityInertialFrame(0);
    aState(4) = initialVelocityInertialFrame(1);
    aState(5) = initialVelocityInertialFrame(2);
    aState(6) = 227;  // Mass [kg] from literature study

    StateAndTime currentStateAndTime(aState);        // Creating the current state class using the namespace and class directly

//////////////////////////////////////////////////////////////////////////////////
/*////////////////////// Testing the Taylor series integrator //////////////////////
//////////////////////////////////////////////////////////////////////////////////

/// Setting the data collection file for TSI and inserting the first values ///

    // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
    const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/01.integrationResults/TSI/";


    // Set output format for matrix output.
    Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

    // Set absolute path to file containing the Taylor Series Coefficients.
    std::string dataAbsolutePath = outputDirectory + "test2TSIstateAndTime.csv";

    // Create a row vector for the storing of the data
    Eigen::MatrixXd outputVector = Eigen::MatrixXd::Zero(1,8); // Create a row vector for the storing of the data

    // Getting the initial conditions for storage
    const tudat::basic_mathematics::Vector7d initialState = currentStateAndTime.getCurrentState();

    // Filling the output vector
    outputVector(0,0) = currentStateAndTime.getCurrentTime();   // Storing the initial time
    outputVector(0,1) = initialState(0);   // Storing the initial x position
    outputVector(0,2) = initialState(1);   // Storing the initial y position
    outputVector(0,3) = initialState(2);   // Storing the initial z position
    outputVector(0,4) = initialState(3);   // Storing the initial x velocity
    outputVector(0,5) = initialState(4);   // Storing the initial y velocity
    outputVector(0,6) = initialState(5);   // Storing the initial z velocity
    outputVector(0,7) = initialState(6);   // Storing the initial MAV mass

    // Storing the data

    std::ifstream ifile(dataAbsolutePath.c_str()); // Check it as an input file

    bool fexists = false;   // Set the default to "It does not exist"

    if (ifile){         // Attempt to open the file


       fexists = true;      // If the file can be opened it must exist

       ifile.close();   // Close the file

    }


    // If so: error and create temporary file, if not: create new file and put data in

    if (fexists == true){


        /// Get time ///

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

//        std::cout<<"currentSec = "<<currentSec<<std::endl;



//        std::string currentRealTime = timeinfo->tm_hour + ";" + timeinfo->tm_min + ":" + timeinfo->tm_sec;

//        std::cout<<"currentRealTime = "<<currentRealTime<<std::endl;

        /// Get time end ///

//        const double currentComputerTime = 11;   // The current CPU time used to make a new file if file already exists (code from stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows)

//        std::string ComputerTimeString;      // Creating a string to hold the CPU time

//        ostringstream convert;          // Stream used for the conversion (as by www.cplusplus.com/articles/D9j2Nwbp/ )

//        convert << currentComputerTime;      // Input the number to be converted

//        ComputerTimeString = convert.str();  // Convert to string and store

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

        std::string ComputerTimeString = currentHour + ":" + currentMin + ":" + currentSec;  // Convert to string and store

        std::string newFileName = "backupTSIFileAtTime_" + ComputerTimeString + ".csv";

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



/// Defining the order and initializing the StepSize class ///

    const int maxOrder = 20;

        StepSize stepSize; // Initializing the stepSize class. THIS SHOULD BE DONE BEFORE THE START OF THE INTEGRATION!!!!!

//        const double currentStepSize = stepSize.getCurrentStepSize();

//        std::cout<<"The current stepSize = "<<currentStepSize<<std::endl;

/// Performing the actual TSI integration ///



        Eigen::VectorXd updatedStateAndTimeVector = performTaylorSeriesIntegrationStep(Mars, MAV, currentStateAndTime, stepSize, maxOrder);
        // This function has the output: updated position, updated velocity, updated mass and updated time

        std::cout<<"updatedStateAndTimeVector = "<<updatedStateAndTimeVector<<std::endl;


        // Check to see if the class has been updated from within the TSI function
//        const double nextStepSize = stepSize.getCurrentStepSize();

//        std::cout<<"The next stepSize = "<<nextStepSize<<std::endl;

/// Storing the values ///

        outputVector = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.

        // Filling the output vector
        outputVector(0,0) = updatedStateAndTimeVector(7);   // Storing the updated time
        outputVector(0,1) = updatedStateAndTimeVector(0);   // Storing the updated x position
        outputVector(0,2) = updatedStateAndTimeVector(1);   // Storing the updated y position
        outputVector(0,3) = updatedStateAndTimeVector(2);   // Storing the updated z position
        outputVector(0,4) = updatedStateAndTimeVector(3);   // Storing the updated x velocity
        outputVector(0,5) = updatedStateAndTimeVector(4);   // Storing the updated y velocity
        outputVector(0,6) = updatedStateAndTimeVector(5);   // Storing the updated z velocity
        outputVector(0,7) = updatedStateAndTimeVector(6);   // Storing the updated MAV mass


        // Check if the file already exists.


        std::ifstream ifile2(dataAbsolutePath.c_str()); // Check it as an input file

        fexists = false;   // Set the default to "It does not exist"

        if (ifile2){         // Attempt to open the file


           fexists = true;      // If the file can be opened it must exist

           ifile2.close();   // Close the file

        }


        // If so: append, if not: create new file and put data in

        if (fexists == true){

            // Export the Taylor Series Coefficients matrix.
            std::ofstream exportFile1;                          // Define the file as an output file


            exportFile1.open(dataAbsolutePath.c_str(),std::ios_base::app);      // Open the file in append mode

            exportFile1 << "\n";                                            // Make sure the new matrix start on a new line

            exportFile1 << outputVector.format( csvFormat ); // Add the new values

            std::cout<<"The file called "<<dataAbsolutePath<<" has been appended"<<std::endl;


            exportFile1.close( );   // Close the file
}
            else{

            std::cerr<<"Error: values could not be stored because storage file does not exist"<<std::endl;
        };

//*/





    return 0;
}


