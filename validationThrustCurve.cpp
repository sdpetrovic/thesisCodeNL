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
 *      160923    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

// This is a file to convert the thrust data to thrust angles

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <numeric> // For vector stuff

#include <quadmath.h>   // Quad precision...

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

// Used for the RKF integrator
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
//#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
//#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
//#include <boost/make_shared.hpp>
//#include <boost/shared_ptr.hpp>
#include <Tudat/Mathematics/NumericalIntegrators/euler.h>
//#include <boost/bind.hpp>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h>

//#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>
#include <Tudat/Mathematics/BasicMathematics/coordinateConversions.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h>
//#include <tudatApplications/thesisProject/referenceFrameTransformationsUpdated.h>

/// Testing the celestial body class ///


#include <thesisProject/celestialBody.h>            // Final version

/// Testing the vehicle class ///

#include <thesisProject/MarsAscentVehicle.h>    // Final version

/// Testing the current state and time and its updater ///

#include <thesisProject/stateAndTime.h>             // Final version


// TSI
/// Testing the auxiliary equations ///
//#include <thesisProject/Auxiliary.h>                // Original test file
#include <thesisProject/AuxiliarySpherical.h>                // Spherical test file


/// Testing the basic recurrence relations ///
#include <thesisProject/projectLibraries/basicRecurrenceRelations.h>               // Original test file

/// Testing all recurrence relations ///
//#include <thesisProject/projectLibraries/allRecurrenceRelations.h>          // Original test file
#include <thesisProject/projectLibraries/allRecurrenceRelationsSpherical.h>          // Spherical test file

/// Testing the stepSize class ///
#include <thesisProject/StepSize.h>             // Original test file

/// Testing the actual Taylor Series integration fucntion ///
//#include <thesisProject/projectLibraries/TaylorSeriesIntegration.h>             // Original test file
#include <thesisProject/projectLibraries/TaylorSeriesIntegrationSpherical.h>             // Original test file

/// Testing the other required functions ///
#include <thesisProject/projectLibraries/otherRequiredFunctions.h>              // Original test file

// State Derivative Function

/// Testing the airTemperature functions ///
#include <thesisProject/projectLibraries/airTemperature.h>      // Original test file

/// Testing the airDensity function ///
#include <thesisProject/projectLibraries/airDensity.h>          // Original test file

/// Testing the ascentDragForce function ///
#include <thesisProject/projectLibraries/ascentDragForce.h>     // Original test file

/// Testing the dragCoefficient function ///
#include <thesisProject/projectLibraries/dragCoefficient.h>     // Original test file

/// Testing the ascentStateDerivativeFunction ///
//#include <thesisProject/projectLibraries/ascentStateDerivativeFunction.h>   // Original test file

/// Testing the ascentStateDerivativeFunctionClass ///
#include <thesisProject/ascentStateDerivativeFunctionClass.h>       // Adapted file from thesisProject/projectLibraries/ascentStateDerivativeFunction.h

// Circularisation
#include <tudat/Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>


// testing
#include <thesisProject/projectLibraries/trajectoryIntegration.h>  // Test the trajectory integration file



int main()

{
std::cout<<setprecision(15)<<"Setting output precision to 15"<<std::endl;
    /// Input file ///

    std::string nameOfFile = "trajectoryOutputDataCase7_3_2016_v33_thrust2.txt"; // Output file from Joel




/// Read the input file ///

    // This requires the complete path in order to work!
    std::string pathToWorkingFile = "/home/stachap/Documents/Thesis/11.ValidationData/";

    std::string completePathToInput = pathToWorkingFile+nameOfFile;

//    std::cout<<"completePathToInput = "<<completePathToInput<<std::endl;

    // This requires the complete path in order to work!
//    std::ifstream inputFile("/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/07.InputFiles/test.txt");
    std::ifstream inputFile(completePathToInput);


    // Read the first line
    std::string firstLine;
    std::getline(inputFile,firstLine);
    std::cout<<firstLine<<std::endl;

    double number1, number2, number3, number4, number5, number6, number7, number8, number9, number10, number11, number12, number13, number14;

    Eigen::MatrixXd inputMatrixValues = Eigen::MatrixXd::Zero(1,14); // Create the first Matrix of the input matrix with all the values


    if (inputFile.is_open())
    {

        int line = 0; // Line to print in vector
        while(inputFile >> number1 >> number2 >> number3 >> number4 >> number5 >> number6 >> number7 >> number8 >> number9 >> number10 >> number11 >> number12 >> number13 >> number14){
//            std::cout<<variable<<" "<<equalSign<<" "<<number<<std::endl;

            inputMatrixValues.conservativeResize(line+1,14); // This resizes the vector and makes it bigger to include all the values in the input file

            inputMatrixValues(line,0) = number1;
            inputMatrixValues(line,1) = number2/1000;
            inputMatrixValues(line,2) = number3/1000;
            inputMatrixValues(line,3) = number4/1000;
            inputMatrixValues(line,4) = number5/1000;
            inputMatrixValues(line,5) = number6/1000;
            inputMatrixValues(line,6) = number7/1000;
            inputMatrixValues(line,7) = number8;
            inputMatrixValues(line,8) = number9;
            inputMatrixValues(line,9) = number10/1000;
            inputMatrixValues(line,10) = deg2rad(number11);  // aoa
            inputMatrixValues(line,11) = number12;
            inputMatrixValues(line,12) = number13;
            inputMatrixValues(line,13) = number14;



            line++;
        }



//        std::cout<<"inputVectorValues"<<'\n'<<inputVectorValues<<std::endl;

        inputFile.close();
    }
    else {std::cout<<"It didn't open..."<<std::endl;}

/// Set input values ///

    std::cout<<"inputMatrixValues = "<<'\n'<<inputMatrixValues<<std::endl;


    celestialBody Mars;



//    const double rotationalVelocityMars = Mars.rotationalVelocity();
    const double primeMeridianAngle = Mars.primeMeridianAngle();
    const double inertialFrameTime = Mars.inertialFrameTime();

//    std::cout<<"rotationalVelocityMars = "<<rotationalVelocityMars<<std::endl;

    Eigen::MatrixXd thrustAngleValuesMatrix = Eigen::MatrixXd::Zero(1,4); // Define an output Matrix

for (int i = 0;i < inputMatrixValues.rows();i++){

    /// Compute current spherical state ///

//        tudat::basic_mathematics::Vector7d currentState = stateAndTime.getCurrentState();

        const double xPosition = inputMatrixValues(i,1);            // x position coordinate definition
        const double yPosition = inputMatrixValues(i,2);            // y position coordinate definition
        const double zPosition = inputMatrixValues(i,3);            // z position coordinate definition
        const double xVelocity = inputMatrixValues(i,4);            // x velocity coordinate definition
        const double yVelocity = inputMatrixValues(i,5);            // y velocity coordinate definition
        const double zVelocity = inputMatrixValues(i,6);            // z velocity coordinate definition
        const double currentTime = inputMatrixValues(i,0);   // current time definition

    // Computations

        const double Radius = sqrt(xPosition*xPosition+yPosition*yPosition+zPosition*zPosition);         // r [km]

        const double currentAltitude = Radius - Mars.bodyReferenceRadius()+0.2; // h [km MOLA]

        const double inertialVelocity = sqrt(xVelocity*xVelocity+yVelocity*yVelocity+zVelocity*zVelocity);       // V_I [km/s]


        const double Latitude = asin(zPosition/Radius);              // delta [rad]

//        std::cout<<"JoelsRotationalVelocityMars = "<<xVelocity/Radius<<std::endl;
        const double rotationalVelocityMars = -xVelocity/Radius; // Set the rotational velocity of mars equal to the one used by Joel


        // Avoid cosine rounding errors
        double cx10;

        if (abs(cos(rotationalVelocityMars*(inertialFrameTime+currentTime)+primeMeridianAngle))<6.2e-17){
            cx10 = 0;
        }
        else {
            cx10 = cos(rotationalVelocityMars*(inertialFrameTime+currentTime)+primeMeridianAngle);
        }



        // Same for sine
        double sx10;

        if (abs(sin(rotationalVelocityMars*(inertialFrameTime+currentTime)+primeMeridianAngle))<6.2e-17){
            sx10 = 0;
        }
        else {
            sx10 = sin(rotationalVelocityMars*(inertialFrameTime+currentTime)+primeMeridianAngle);
        }



        const double rotationalXposition = cx10*xPosition+sx10*yPosition;       // x_R [km]

        const double rotationalYposition = -sx10*xPosition+cx10*yPosition;      // y_R [km]

        const double rotationalLongitude = atan2(rotationalYposition,rotationalXposition); // tau [rad]
        /// Debug ///
        //std::cout<<"Here 2"<<std::endl;
        /// Debug ///


//        double rotationalVelocity_;   // V_G [km/s] (placeholder)
//        // Deal with rounding errors
//        if (inertialVelocity*inertialVelocity+rotationalVelocityMars*rotationalVelocityMars*(xPosition*xPosition+yPosition*yPosition)+2.0*rotationalVelocityMars*(xVelocity*yPosition-yVelocity*xPosition) <= 0.0){
//            rotationalVelocity_ = 0.0;
//        }
//        else {
//         rotationalVelocity_ = sqrt(inertialVelocity*inertialVelocity+rotationalVelocityMars*rotationalVelocityMars*(xPosition*xPosition+yPosition*yPosition)+2.0*rotationalVelocityMars*(xVelocity*yPosition-yVelocity*xPosition)); // V_G [km/s]
//        }

//        const double rotationalVelocity = rotationalVelocity_;      // V_G [km/s] (actual parameter)
        const double rotationalVelocity = inputMatrixValues(i,9); // V_G [km/s] taken from the provided data sheet

         /// Debug ///
//        std::cout<<"V_I^2*+Omega_M^2*(x1^2+x2^2)+2.0*Omega_M*(x4*x2-x5*x1) = "<<inertialVelocity*inertialVelocity+rotationalVelocityMars*rotationalVelocityMars*(xPosition*xPosition+yPosition*yPosition)+2.0*rotationalVelocityMars*(xVelocity*yPosition-yVelocity*xPosition)<<std::endl;
        //std::cout<<"Here 3"<<std::endl;
        /// Debug ///

        // Avoid cosine rounding errors
        double cx11;

        if (abs(cos(rotationalLongitude))<6.2e-17){
            cx11 = 0;
        }
        else {
            cx11 = cos(rotationalLongitude);
        }


        double cx12;

        if (abs(cos(Latitude))<6.2e-17){
            cx12 = 0;
        }
        else {
            cx12 = cos(Latitude);
        }

        // Same for sine
        double sx11;

        if (abs(sin(rotationalLongitude))<6.2e-17){
            sx11 = 0;
        }
        else {
            sx11 = sin(rotationalLongitude);
        }


        double sx12;

        if (abs(sin(Latitude))<6.2e-17){
            sx12 = 0;
        }
        else {
            sx12 = sin(Latitude);
        }




        const double verticalXvelocity = (xVelocity+rotationalVelocityMars*yPosition)*(sx12*sx10*sx11-
                                                                                       cx10*cx11*sx12)+
                (yVelocity-rotationalVelocityMars*xPosition)*(-cx11*sx12*sx10-
                                                              cx10*sx12*sx11)+zVelocity*cx12;   // Vx_V [km/s]

        const double verticalYvelocity = (xVelocity+rotationalVelocityMars*yPosition)*(-cx11*sx10-
                                                                                       cx10*sx11)+
                (yVelocity-rotationalVelocityMars*xPosition)*(cx10*cx11-
                                                              sx10*sx11);     // Vy_V [km/s]

        const double verticalZvelocity = (xVelocity+rotationalVelocityMars*yPosition)*(cx12*sx10*sx11-
                                                                                       cx10*cx11*cx12)+
                (yVelocity-rotationalVelocityMars*xPosition)*(-cx11*cx12*sx10-
                                                              cx10*cx12*sx11)-zVelocity*sx12;   // Vz_V [km/s]

        double rotationalAzimuth_;
        double rotationalFlightPathAngle_;          // gamma_G [rad] (placeholder)


            rotationalAzimuth_ = atan2(verticalYvelocity,verticalXvelocity);    // chi_G [rad]

            // Avoid singularities
            if (rotationalVelocity == 0.0){
                rotationalFlightPathAngle_ = tudat::mathematical_constants::LONG_PI/2.0;
            }
            else if (verticalZvelocity/rotationalVelocity >= 1.0 || verticalZvelocity/rotationalVelocity-1 >= -1E-15){  // Compensate for rounding errors
                rotationalFlightPathAngle_ = -asin(1.0);
            }
            else if (verticalZvelocity/rotationalVelocity <= -1.0 || verticalZvelocity/rotationalVelocity+1 <= 1E-15){ // Compensate for rounding errors
                rotationalFlightPathAngle_ = -asin(-1.0);
            }
            else {
            rotationalFlightPathAngle_ = -asin(verticalZvelocity/rotationalVelocity);   // gamma_G [rad]

            }

//std::cout<<"Vzv = "<<verticalZvelocity<<std::endl;

        const double rotationalAzimuth = rotationalAzimuth_;        // chi_G [rad]


        const double rotationalFlightPathAngle = rotationalFlightPathAngle_; // gamma_G [rad] (actual parameter)


        /// Check output ///

//        std::cout<<"Radius = "<<Radius<<std::endl;
//        std::cout<<"rotationalLongitude = "<<rotationalLongitude<<std::endl;
//        std::cout<<"Latitude = "<<Latitude<<std::endl;
//        std::cout<<"rotationalVelocity = "<<rotationalVelocity<<std::endl;
//        std::cout<<"rotationalAzimuth = "<<rotationalAzimuth<<std::endl;
//        std::cout<<"rotationalFlightPathAngle = "<<rotationalFlightPathAngle<<std::endl;


/////////////////////////////////////////////////////////////////////// Transformations ///////////////////////////////////////////////////////////////////////


////    std::cout<<"inertialVelocityR"<<inertialVelocityR<<std::endl;


    /// Thrust directions in the inertial frame ///

    Eigen::Vector3d inertialThrustDirection;

    inertialThrustDirection(0) = inputMatrixValues(i,11); // x-unit vector
    inertialThrustDirection(1) = inputMatrixValues(i,12); // y-unit vector
    inertialThrustDirection(2) = inputMatrixValues(i,13); // z-unit vector

//    std::cout<<"inertialThrustDirection = "<<inertialThrustDirection<<std::endl;

    /// Thrust direction in the rotational frame ///
    const Eigen::Vector3d rotationalThrustDirection = tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix(rotationalVelocityMars*currentTime)*inertialThrustDirection;

    /// Thrust direction in the vertical frame ///
    const Eigen::Vector3d verticalThrustDirection = tudat::reference_frames::getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(rotationalLongitude, Latitude)*rotationalThrustDirection;


    /// Thrust direction in the aerodynamic frame ///
    const Eigen::Vector3d aerodynamicThrustDirection = tudat::reference_frames::getLocalVerticalFrameToTrajectoryTransformationMatrix(rotationalFlightPathAngle, rotationalAzimuth)*verticalThrustDirection;


    /// Thrust direction in the body frame ///
    const Eigen::Vector3d bodyThrustDirection = tudat::reference_frames::getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(inputMatrixValues(i,10),0.0)*aerodynamicThrustDirection;

//    std::cout<<"bodyThrustDirection = "<<bodyThrustDirection<<std::endl;

/////////////////////////////////////////////////////////////////////// Compute the angles ///////////////////////////////////////////////////////////////////////

    const double thrustAzimuthAngle = atan2(bodyThrustDirection(1),bodyThrustDirection(0));
    const double thrustElevationAngle = atan2(-bodyThrustDirection(2),sqrt(bodyThrustDirection(0)*bodyThrustDirection(0)+bodyThrustDirection(1)*bodyThrustDirection(1)));

    std::cout<<"thrustAzimuthAngle [deg] = "<<rad2deg(thrustAzimuthAngle)<<std::endl;
    std::cout<<"thrustElevationAngle [deg] = "<<rad2deg(thrustElevationAngle)<<std::endl;
    std::cout<<"combinedAngle [deg] = "<<rad2deg(acos(cos(thrustAzimuthAngle)*cos(thrustElevationAngle)))<<std::endl;
    std::cout<<"angleOfAttack = "<<rad2deg(inputMatrixValues(i,10))<<std::endl;

    // Check to see if the rotation will actually produce a thrust in just the x-direction in the propulsion frame
    const Eigen::Vector3d propulsionThrustDirection = getBodyToPropulsionFrameTransformationQuaternion(thrustAzimuthAngle, thrustElevationAngle)*bodyThrustDirection;

//    std::cout<<"propulsionThrustDirection = "<<propulsionThrustDirection<<std::endl;


    thrustAngleValuesMatrix.conservativeResize(i+1,4); // Resize the matrix to fit all the values

    thrustAngleValuesMatrix(i,0) = currentTime;
    thrustAngleValuesMatrix(i,1) = currentAltitude;
    thrustAngleValuesMatrix(i,2) = rad2deg(thrustAzimuthAngle);
    thrustAngleValuesMatrix(i,3) = rad2deg(thrustElevationAngle);

/////////////////////////////////////////////////////////////////////// Take two ///////////////////////////////////////////////////////////////////////

    /// Inertial velocity in the rotating frame ///

            const double inertialXvelocityR = xVelocity+rotationalVelocityMars*yPosition;
            const double inertialYvelocityR = yVelocity-rotationalVelocityMars*xPosition;

        Eigen::Vector3d inertialVelocityR;

        inertialVelocityR(0) = inertialXvelocityR;
        inertialVelocityR(1) = inertialYvelocityR;
        inertialVelocityR(2) = zVelocity;

        Eigen::Vector3d rotationalVelocityVector = tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix(rotationalVelocityMars*currentTime)*inertialVelocityR;


//        std::cout<<"rotationalVelocityVector = "<<rotationalVelocityVector<<std::endl;

        double innerProductVectors = 0;

        for (int j=0;j<3;j++){
            innerProductVectors+=inertialThrustDirection(j)*rotationalVelocityVector(j);
        }
//std::cout<<"rotationalVelocity = "<<rotationalVelocity<<std::endl;
//std::cout<<"innerProductVectors = "<<innerProductVectors<<std::endl;

//        const double thrustAngle = rad2deg(acos((inertialThrustDirection*rotationalVelocityVector)/rotationalVelocity));
        const double thrustAngle = rad2deg(acos(innerProductVectors/rotationalVelocity));

        std::cout<<"thrustAngle [deg] = "<<thrustAngle<<std::endl;


} // end of for loop






/// Save the angles ///

//    /*        // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
        const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/04.VerificationAndValidation/";


        // Set output format for matrix output.
        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );


//            /// Get time ///

//            time_t rawtime;
//            struct tm * timeinfo;

//            time ( &rawtime );
//            timeinfo = localtime ( &rawtime );
//    //        printf ( "Current local time and date: %s", asctime (timeinfo) );         // (from stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c and from stackoverflow.com/questions/1442116/how-to-get-date-and-time-value-in-c-program)
//    //        std::cout<<"time = "<<time ( &rawtime )<<std::endl;
//    //        std::cout<< (timeinfo->tm_year+1900) << "-"
//    //                 << (timeinfo->tm_mon + 1) << "-"
//    //                 << timeinfo->tm_mday << " "
//    //                 << timeinfo->tm_hour << ":"
//    //                 << timeinfo->tm_min << ":"
//    //                 << timeinfo->tm_sec <<std::endl;

//            ostringstream ConvertYear;
//            ConvertYear << (timeinfo->tm_year+1900);
//            std::string currentYear = ConvertYear.str();

//            ostringstream ConvertMonth;
//            ConvertMonth << (timeinfo->tm_mon+1);
//            std::string currentMonth = ConvertMonth.str();

//            ostringstream ConvertDay;
//            ConvertDay << timeinfo->tm_mday;
//            std::string currentDay = ConvertDay.str();

//            ostringstream ConvertHour;
//            ConvertHour << timeinfo->tm_hour;
//            std::string currentHour = ConvertHour.str();

//    //        std::cout<<"currentHour = "<<currentHour<<std::endl;

//            ostringstream ConvertMin;
//            ConvertMin << timeinfo->tm_min;
//            std::string currentMin = ConvertMin.str();

//    //        std::cout<<"currentMin = "<<currentMin<<std::endl;

//            ostringstream ConvertSec;
//            ConvertSec << timeinfo->tm_sec;
//            std::string currentSec = ConvertSec.str();

//    //        std::cout<<"currentSec = "<<currentSec<<std::endl;

//            // Making sure each one of the representation has at two numbers
//            if(currentMonth.size() == 1){
//                currentMonth = "0" + currentMonth;
//            }


//            if(currentDay.size() == 1){
//                currentDay = "0" + currentDay;
//            }


//            if(currentSec.size() == 1){
//                currentSec = "0" + currentSec;
//            }

//            if (currentMin.size() == 1){
//                currentMin = "0" + currentMin;
//            }

//            if (currentHour.size() == 1){
//                currentHour = "0" + currentHour;
//            }


//    //        std::cout<<"The length of currentSec = "<<currentSec.size()<<" and the value = "<<currentSec<<std::endl;

//            std::string ComputerTimeString = currentYear + "-" + currentMonth + "-" + currentDay + "_" + currentHour + ":" + currentMin + ":" + currentSec;  // Convert to string and store

//            std::string newFileName = "SphericalTaylorSeriesCoefficientsFileAtDateAndTime_" + ComputerTimeString + ".csv";

        std::string newFileName = "ThrustAngles.csv";

//    std::cerr<<"The file name that you have chosen already exists, a new file with name "<<newFileName<<" will be created to store the data for now"<<std::endl;

//    // Set new absolute path to file containing the data.
//    dataAbsolutePathTSI = outputDirectoryTSI + newFileName;

    // Set absolute path to file containing the Taylor Series Coefficients.
    const std::string thrustAnglesAbsolutePath = outputDirectory + newFileName;


        // Set absolute path to file containing the Taylor Series Coefficients.
//        const std::string taylorSeriesCoefficientsAbsolutePath = outputDirectory + "test1TaylorSeriesCoefficients(bugSearch31-05-2016).csv";



        // Check if the file already exists.


        std::ifstream ifile(thrustAnglesAbsolutePath.c_str()); // Check it as an input file

        bool fexists = false;   // Set the default to "It does not exist"

        if (ifile){         // Attempt to open the file


           fexists = true;      // If the file can be opened it must exist

           ifile.close();   // Close the file

        }


        // If so: append, if not: create new file and put data in

        if (fexists == true){

//            // Export the Taylor Series Coefficients matrix.
//            std::ofstream exportFile1;                          // Define the file as an output file


//            exportFile1.open(thrustAnglesAbsolutePath.c_str(),std::ios_base::app);      // Open the file in append mode

//            exportFile1 << "\n";                                            // Make sure the new matrix start on a new line

//            exportFile1 << thrustAngleValuesMatrix.format( csvFormat ); // Add the new values

////            std::cout<<"The file called "<<taylorSeriesCoefficientsAbsolutePath<<" has been appended"<<std::endl;


//            exportFile1.close( );   // Close the file

            std::cout<<"The current thrust angle file has to be deleted first!"<<std::endl;
}
            else{

            // Export the Taylor Series Coefficients thrustAnglesAbsolutePath.
            std::ofstream exportFile1( thrustAnglesAbsolutePath.c_str( ) ); // Make the new file
            std::cout<<"New file called "<<thrustAnglesAbsolutePath<<" has been created"<<std::endl;
            exportFile1 << thrustAngleValuesMatrix.format( csvFormat );          // Store the new values
            exportFile1.close( );   // Close the file
        };

        //*/

    return 0;
}


