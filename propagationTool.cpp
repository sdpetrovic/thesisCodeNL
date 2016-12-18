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
 *      160919    S.D. Petrovic     File created
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

 /*   /// Running instructions ///

    Whenever a new run is made, make sure these things are correct:
    1. - The case that is being run
    2. - The currentVariable string is correct
    3. - The current variable is deactivated in the /// Set input values /// section
    4. - The desired end time is chosen
    5. - The planet is rotating or not
    6. - The correct loop is set for the current variable
    7. - The current variable is printed in the first column of the output matrix (see section /// Storing the output TSI)

  //*/  /// Running instructions ///

    /// Input file ///

    int Case = 2;  // Case 1 = Ryan, Case 2 = Joel

    std::string nameOfFile; // Declare the name of the file string
    std::string caseName; // Declare the case name used for the file location
    const std::string currentVariable = "headingAngle"; // Delcaare the current variable that is changed

    /* List of possible variable names
     *
     *  - Order
     *  - ErrorTolerance
     *  - launchAltitude
     *  - launchLatitude
     *  - launchLongitude
     *  - Position         (Such as (0,0) (0,90) and (90,0) latitude and longitude)
     *  - FPA
     *  - headingAngle
     *  - groundVelocity
     *  - normalRun
     *
     *
     *
     *
     *
     *
     * //*/

//    std::string nameOfFile = "test.txt"; // Test file

    if (Case == 1){
        nameOfFile = "test1Ryan.txt"; // Test file from Ryan
        caseName = "Case1"; // Test folder used for output associated with this case
    }
    else if (Case == 2){
        nameOfFile = "test2Joel.txt"; // Test file from Joel
        caseName = "Case2"; // Test folder used for output associated with this case

    }






/// Read the input file ///

    // This requires the complete path in order to work!
    std::string pathToWorkingFile = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/07.InputFiles/";

    std::string completePathToInput = pathToWorkingFile+nameOfFile;

//    std::cout<<"completePathToInput = "<<completePathToInput<<std::endl;

    // This requires the complete path in order to work!
//    std::ifstream inputFile("/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/07.InputFiles/test.txt");
    std::ifstream inputFile(completePathToInput);


    // Read the first line
    std::string firstLine, variable, equalSign;
    std::getline(inputFile,firstLine);
    std::cout<<firstLine<<std::endl;

    double number;

    Eigen::VectorXd inputVectorValues = Eigen::VectorXd::Zero(1); // Create the first vector of the input matrix with all the values


    if (inputFile.is_open())
    {

        int line = 0; // Line to print in vector
        while(inputFile >> variable >> equalSign >> number){
//            std::cout<<variable<<" "<<equalSign<<" "<<number<<std::endl;

            inputVectorValues.conservativeResize(line+1); // This resizes the vector and makes it bigger to include all the values in the input file

            inputVectorValues(line) = number;

            line++;
        }



//        std::cout<<"inputVectorValues"<<'\n'<<inputVectorValues<<std::endl;

        inputFile.close();
    }
    else {std::cout<<"It didn't open..."<<std::endl;}

/// Set input values ///

//    const double desiredOrbitalAltitude =   inputVectorValues(0);
    const double desiredOrbitalAltitude = 20000; // Case for determination of when zero FPA is reached

    const double desiredInclinationDeg =    inputVectorValues(1);
    const double initialAltitude =          inputVectorValues(2);
    const double initialLatitudeDeg =       inputVectorValues(3);
    const double initialLongitudeDeg =      inputVectorValues(4);

    const double FlightPathAngleDeg =       inputVectorValues(5);
//    const double FlightPathAngleDeg = 89;

//    const double HeadingAngleDeg =          inputVectorValues(6);

    const double initialGroundVelocity =    inputVectorValues(7);
//    const double initialGroundVelocity = 1e-6;

    const double massMAV =                  inputVectorValues(8);
    const double thrust =                   inputVectorValues(9);
    const double specificImpulse =          inputVectorValues(10);
    const double initialBurnTime =          inputVectorValues(11);
    const double constantThrustElevationAngle = inputVectorValues(12);
    const double constantThrustAzimuthAngle = inputVectorValues(13);

    const int maxOrder =                    inputVectorValues(14);
//    const int maxOrder = 20;

    const double chosenLocalErrorTolerance = inputVectorValues(15);
//    const double chosenLocalErrorTolerance = 1e-15;

    const double chosenStepSize =           inputVectorValues(16);


//    const double setEndTime =               inputVectorValues(17);

//    const double setEndTime = 20000;        // Case for determination of when zero FPA is reached


    const double setEndTime = 1796;         // End time for tests for Ryans case 1 (For state accuracy)
//    const double setEndTime = 789;          // End time for tests for Ryans case 1 (For state accuracy) for non rotating planet. Also used for comparison between rotating and non-rotating
//    const double setEndTime = 1543;
//    const double setEndTime = 1795;
//    const double setEndTime = 1436;
//    const double setEndTime = 189;           // End time to determine error in order


//    const double setEndTime = 876;        // End time for tests for Joels case 2. Also used for comparison between rotating and non-rotating
//    const double setEndTime = 783;          // For Latitude
//    const double setEndTime = 856;          // For altitude
//    const double setEndTime = 873;        // For Longitude
//    const double setEndTime = 928;          // Time for non rotating planet. Latitude and longitude


    const double RKFinitiaterTime =         inputVectorValues(18);

    const bool rotatingPlanet =             inputVectorValues(19);
//    const bool rotatingPlanet = false;

    const bool GravityAcc =                 inputVectorValues(20);
    const bool ThrustAcc =                  inputVectorValues(21);
    const bool DragAcc =                    inputVectorValues(22);
    const bool comparison =                 inputVectorValues(23);


    std::string rotOrNot; // Define the folder, either rotating Mars or not

    if (rotatingPlanet == true){
        rotOrNot = "rotatingMars";
    }
    else if (rotatingPlanet == false){
        rotOrNot = "nonRotatingMars";
    }


    //////////// Variables for storing of the data /////////////////

 /*  //// Order version 1 ////
//    const std::string currentVariable = "Order";
    const int maximumFinalOrder = 29;
    const int initialOrder = 21;
    const int orderStepSize = 2;
//    const double numberOfRuns = (maximumFinalOrder+1-initialOrder)/orderStepSize; // Number of runs that order will be changed and evaluated
    const double numberOfRuns = (maximumFinalOrder-initialOrder)/orderStepSize+1;
    Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
    Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
    int run = 0;    // Initual run counter to be used in the storing of the output


    /// For loop ///

    for (int maxOrder = initialOrder; maxOrder<maximumFinalOrder+1; maxOrder = maxOrder+orderStepSize){
        std::cout<<"Max order = "<<maxOrder<<std::endl;
//*/    //// Order version 1////


  /*      //// Order version 2 ////
     //    const std::string currentVariable = "Order";
//         const int maximumFinalOrder = 32;
//         const int initialOrder = 5;
//         const int orderStepSize = 1;
//         const double numberOfRuns = (maximumFinalOrder+1-initialOrder)/orderStepSize; // Number of runs that order will be changed and evaluated


        const double startOrder = 5.0;
        const double endOrder = 31.0;    // 100 for Case 1 and 31 for Case 2

        const int removeNumbers = 0;   // REMOVE THE -3 AFTER TESTING FOR 25, 27 AND 29!

         const double numberOfRuns = endOrder-startOrder+1-removeNumbers;
         std::cout<<"Number of runs = "<<numberOfRuns<<std::endl;

         Eigen::VectorXd maxOrderVector = Eigen::VectorXd::Zero(numberOfRuns);

//         maxOrderVector(0) = 3;
//         maxOrderVector(1) = 4;
//         maxOrderVector(2) = 5;
//         maxOrderVector(3) = 6;
//         maxOrderVector(4) = 7;
//         maxOrderVector(5) = 8;
//         maxOrderVector(6) = 9;
//         maxOrderVector(7) = 10;
//         maxOrderVector(8) = 11;
//         maxOrderVector(9) = 12;
//         maxOrderVector(10) = 13;
//         maxOrderVector(11) = 14;
//         maxOrderVector(12) = 15;
//         maxOrderVector(13) = 16;
//         maxOrderVector(14) = 17;
//         maxOrderVector(15) = 18;
//         maxOrderVector(16) = 19;
//         maxOrderVector(17) = 20;
//         maxOrderVector(18) = 21;
//         maxOrderVector(19) = 22;
//         maxOrderVector(20) = 23;
//         maxOrderVector(21) = 24;
//         maxOrderVector(22) = 25;
//         maxOrderVector(23) = 26;
//         maxOrderVector(24) = 27;
//         maxOrderVector(25) = 28;
//         maxOrderVector(26) = 29;
//         maxOrderVector(27) = 30;

//         int kk = startOrder;
//         for (int j = 0; j<numberOfRuns; j++){ // More efficient way to fill the vector

//             if (kk==25 || kk==27 || kk==29){
//                kk++;
//                maxOrderVector(j) = kk;
//                kk++;
//             }
//             else {
//             maxOrderVector(j) = kk;
//             kk++;
//             }
//         }


         int kk = startOrder;
         for (int j = 0; j<numberOfRuns; j++){ // More efficient way to fill the vector

             maxOrderVector(j) = kk;
             kk++;

         }







         double maxOrder = 0; // The actual parameter

         Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
         Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
         int run = 0;    // Initual run counter to be used in the storing of the output


         /// For loop ///

         for (int i = 0; i<numberOfRuns; i++){
             maxOrder = maxOrderVector(i);
             std::cout<<"Max order = "<<maxOrder<<std::endl;
     //*/    //// Order version 2////

  /*   //// Error tolerance ////
//        const std::string currentVariable = "ErrorTolerance";
//        const double maximumFinalTolerance = 1e-15;
//        const double initialTolerance = 1e-5;
    const double numberOfRuns = 11; // Number of runs that order will be changed and evaluated

        Eigen::VectorXd chosenLocalErrorToleranceVector = Eigen::VectorXd::Zero(numberOfRuns);

        chosenLocalErrorToleranceVector(0) = 1e-5;
        chosenLocalErrorToleranceVector(1) = 1e-6;
        chosenLocalErrorToleranceVector(2) = 1e-7;
        chosenLocalErrorToleranceVector(3) = 1e-8;
        chosenLocalErrorToleranceVector(4) = 1e-9;
        chosenLocalErrorToleranceVector(5) = 1e-10;
        chosenLocalErrorToleranceVector(6) = 1e-11;
        chosenLocalErrorToleranceVector(7) = 1e-12;
        chosenLocalErrorToleranceVector(8) = 1e-13;
        chosenLocalErrorToleranceVector(9) = 1e-14;
        chosenLocalErrorToleranceVector(10) = 1e-15;



        double chosenLocalErrorTolerance = 0; // The actual parameter


//        const double orderStepSize = 0.1;

        Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
        Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
        int run = 0;    // Initual run counter to be used in the storing of the output


        /// For loop ///

        for (int i = 0; i<numberOfRuns; i++){
            chosenLocalErrorTolerance = chosenLocalErrorToleranceVector(i);
            std::cout<<"chosenLocalErrorTolerance = "<<chosenLocalErrorTolerance<<std::endl;
    //*/    //// Error Tolerance ////


 /*  //// Initial Altitude ////
//    const double numberOfRuns = 3;


//    // 3 data points
//    const double startAltitude = -0.5;
//    const double endAltitude = 0.5;
//    const double altitudeStepSize = 0.5;

    // complete range
    const double startAltitude = -0.6;
    const double endAltitude = 0.5;
    const double altitudeStepSize = 0.1;

     const double numberOfRuns = (endAltitude-startAltitude)/altitudeStepSize+1;
     std::cout<<"Number of runs = "<<numberOfRuns<<std::endl;

     Eigen::VectorXd initialAltitudeVector = Eigen::VectorXd::Zero(numberOfRuns);
    double initialStep = 0;
         for (int j = 0; j<numberOfRuns; j++){ // More efficient way to fill the vector
             initialAltitudeVector(j) = startAltitude+initialStep;
             initialStep = initialStep + altitudeStepSize;
             std::cout<<"initialStep = "<<initialStep<<std::endl;
         }


    double initialAltitude = 0; // The actual parameter


    Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
    Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
    int run = 0;    // Initual run counter to be used in the storing of the output


    /// For loop ///

    for (int i = 0; i<numberOfRuns; i++){
        initialAltitude = initialAltitudeVector(i);
        std::cout<<"initialAltitude = "<<initialAltitude<<std::endl;




    //*/    //// Initial Altitude ////


/*    //// Initial Position ////
    const double numberOfRuns = 5;

    Eigen::MatrixXd initialPositionMatrix = Eigen::MatrixXd::Zero(numberOfRuns,2);

    initialPositionMatrix(0,0) = 0.0; // initialLatitudeDeg
    initialPositionMatrix(0,1) = 0.0; // initialLongitudeDeg

//        // Case 1
//        initialPositionMatrix(1,0) = 30.0;
//        initialPositionMatrix(1,1) = 0.0;

//        initialPositionMatrix(2,0) = 0.0;
//        initialPositionMatrix(2,1) = 30.0;

//        initialPositionMatrix(3,0) = 30.0;
//        initialPositionMatrix(3,1) = 30.0;

//        initialPositionMatrix(4,0) = 45.0;
//        initialPositionMatrix(4,1) = 45.0;

//        initialPositionMatrix(5,0) = 0.0;
//        initialPositionMatrix(5,1) = 90.0;

      // Case 2
    initialPositionMatrix(1,0) = -30.0;
    initialPositionMatrix(1,1) = 0.0;

    initialPositionMatrix(2,0) = 0.0;
    initialPositionMatrix(2,1) = 30.0;

    initialPositionMatrix(3,0) = -30.0;
    initialPositionMatrix(3,1) = 30.0;

    initialPositionMatrix(4,0) = 0.0;
    initialPositionMatrix(4,1) = 90.0;



    double initialLatitudeDeg = 0; // The actual parameter
    double initialLongitudeDeg = 0; // The actual parameter



    Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
    Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
    int run = 0;    // Initual run counter to be used in the storing of the output


    /// For loop ///

    for (int i = 0; i<numberOfRuns; i++){
        initialLatitudeDeg = initialPositionMatrix(i,0);
        initialLongitudeDeg = initialPositionMatrix(i,1);
        std::cout<<"initialLatitudeDeg = "<<initialLatitudeDeg<<std::endl;
        std::cout<<"initialLongitudeDeg = "<<initialLongitudeDeg<<std::endl;

    //*/  //// Initial Position ///


/*    //// Initial Latitude ////
    const double numberOfRuns = 11;

    Eigen::VectorXd initialLatitudeVector = Eigen::VectorXd::Zero(numberOfRuns);

    // Case 1
    initialLatitudeVector(0) = 0;
    initialLatitudeVector(1) = 10;
    initialLatitudeVector(2) = 20;
    initialLatitudeVector(3) = 30;
    initialLatitudeVector(4) = 40;
    initialLatitudeVector(5) = 50;
    initialLatitudeVector(6) = 60;
    initialLatitudeVector(7) = 70;
    initialLatitudeVector(8) = 80;
    initialLatitudeVector(9) = 89.9;
    initialLatitudeVector(10) = 90;

//    // Case 2
//    initialLatitudeVector(0) = 0;
//    initialLatitudeVector(1) = -10;
//    initialLatitudeVector(2) = -20;
//    initialLatitudeVector(3) = -30;
//    initialLatitudeVector(4) = -40;
//    initialLatitudeVector(5) = -50;
//    initialLatitudeVector(6) = -60;
//    initialLatitudeVector(7) = -70;
//    initialLatitudeVector(8) = -80;
//    initialLatitudeVector(9) = -89.9;
//    initialLatitudeVector(10) = -90;

    double initialLatitudeDeg = 0; // The actual parameter
    double initialLongitudeDeg = 0; // The actual parameter



    Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
    Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
    int run = 0;    // Initual run counter to be used in the storing of the output


    /// For loop ///

    for (int i = 0; i<numberOfRuns; i++){
        initialLatitudeDeg = initialLatitudeVector(i);
        std::cout<<"initialLatitudeDeg = "<<initialLatitudeDeg<<std::endl;
//        std::cout<<"initialLongitudeDeg = "<<initialLongitudeDeg<<std::endl;

    //*/  //// Initial Latitude ///


 /*     //// Initial Longitude ////
        const double numberOfRuns = 12;

        Eigen::VectorXd initialLongitudeVector = Eigen::VectorXd::Zero(numberOfRuns);


        initialLongitudeVector(0) = 0;
        initialLongitudeVector(1) = 30;
        initialLongitudeVector(2) = 60;
        initialLongitudeVector(3) = 90;
        initialLongitudeVector(4) = 120;
        initialLongitudeVector(5) = 150;
        initialLongitudeVector(6) = 180;
        initialLongitudeVector(7) = 210;
        initialLongitudeVector(8) = 240;
        initialLongitudeVector(9) = 270;
        initialLongitudeVector(10) = 300;
        initialLongitudeVector(11) = 330;




//        double initialLatitudeDeg = 0; // The actual parameter
        double initialLongitudeDeg = 0; // The actual parameter



        Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
        Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
        int run = 0;    // Initual run counter to be used in the storing of the output


        /// For loop ///

        for (int i = 0; i<numberOfRuns; i++){
            initialLongitudeDeg = initialLongitudeVector(i);
            std::cout<<"initialLongitudeDeg = "<<initialLongitudeDeg<<std::endl;
    //        std::cout<<"initialLongitudeDeg = "<<initialLongitudeDeg<<std::endl;

        //*/  //// Initial Longitude ///


 /*      //// FPA ////
           const double numberOfRuns = 7;

           Eigen::VectorXd FPAVector = Eigen::VectorXd::Zero(numberOfRuns);

//           // Case 1
//           FPAVector(0) = 89;
//           FPAVector(1) = 88;
//           FPAVector(2) = 87;




           // Case 2
           FPAVector(0) = 89;
           FPAVector(1) = 88;
           FPAVector(2) = 87;
           FPAVector(3) = 86;
           FPAVector(4) = 85;
           FPAVector(5) = 84;
           FPAVector(6) = 83;







           double FlightPathAngleDeg = 0; // The actual parameter



           Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
           Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
           int run = 0;    // Initual run counter to be used in the storing of the output


           /// For loop ///

           for (int i = 0; i<numberOfRuns; i++){
               FlightPathAngleDeg = FPAVector(i);
               std::cout<<"FlightPathAngleDeg = "<<FlightPathAngleDeg<<std::endl;

           //*/  //// FPA ///

    //// Heading angle ////
               const double numberOfRuns = 12;

               Eigen::VectorXd headingAngleVector = Eigen::VectorXd::Zero(numberOfRuns);

//               headingAngleVector(0) = 0;
//               headingAngleVector(1) = 45;
//               headingAngleVector(2) = 90;
//               headingAngleVector(3) = 135;
//               headingAngleVector(4) = 180;


//               headingAngleVector(0) = 0;
//               headingAngleVector(1) = 45;
//               headingAngleVector(2) = 90;
//               headingAngleVector(3) = 135;
//               headingAngleVector(4) = 180;
//               headingAngleVector(5) = 225;
//               headingAngleVector(6) = 270;
//               headingAngleVector(7) = 315;



               headingAngleVector(0) = 0;
               headingAngleVector(1) = 30;
               headingAngleVector(2) = 60;
               headingAngleVector(3) = 90;
               headingAngleVector(4) = 120;
               headingAngleVector(5) = 150;
               headingAngleVector(6) = 180;
               headingAngleVector(7) = 210;
               headingAngleVector(8) = 240;
               headingAngleVector(9) = 270;
               headingAngleVector(10) = 300;
               headingAngleVector(11) = 330;




               double HeadingAngleDeg = 0; // The actual parameter



               Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
               Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
               int run = 0;    // Initual run counter to be used in the storing of the output


               /// For loop ///

               for (int i = 0; i<numberOfRuns; i++){
                   HeadingAngleDeg = headingAngleVector(i);
                   std::cout<<"HeadingAngleDeg = "<<HeadingAngleDeg<<std::endl;

               //*/  //// Heading angle ///

/*    //// Ground velocity ////
                   const double numberOfRuns = 6;

                   Eigen::VectorXd groundVelocityVector = Eigen::VectorXd::Zero(numberOfRuns);

                   groundVelocityVector(0) = 1e-6;
                   groundVelocityVector(1) = 1e-5;
                   groundVelocityVector(2) = 1e-4;
                   groundVelocityVector(3) = 1e-3;
                   groundVelocityVector(4) = 1e-2;
                   groundVelocityVector(5) = 1e-1;







                   double initialGroundVelocity = 0; // The actual parameter



                   Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
                   Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
                   int run = 0;    // Initual run counter to be used in the storing of the output


                   /// For loop ///

                   for (int i = 0; i<numberOfRuns; i++){
                       initialGroundVelocity = groundVelocityVector(i);
                       std::cout<<"initialGroundVelocity = "<<initialGroundVelocity<<std::endl;

                   //*/  //// Ground velocity ///


/*  //// Normal run ////

    const double numberOfRuns = 1; // These are the actual number of runs that are done
    const double maxOrder = 22;

    Eigen::VectorXd differentEndTimes = Eigen::VectorXd::Zero(numberOfRuns);

    double kk = 17.0; // 14.25
    double stepSizeLoop = 1.0;
    for (int j = 0; j<numberOfRuns; j++){ // More efficient way to fill the vector

        differentEndTimes(j) = kk;
        kk = kk+stepSizeLoop;

    }




    double setEndTime = 0; // The actual parameter

    Eigen::MatrixXd outputMatrixTSI = Eigen::MatrixXd::Zero(numberOfRuns,28); // Creating the output matrix for TSI
    Eigen::MatrixXd outputMatrixRKF = Eigen::MatrixXd::Zero(numberOfRuns,14); // Creating the ouptut matrix for RKF
    int run = 0;    // Initual run counter to be used in the storing of the output

    for (int i = 0; i<numberOfRuns; i++){
        setEndTime = differentEndTimes(i);
        std::cout<<"setEndTime = "<<setEndTime<<std::endl;


 //*/   //// Normal run ////



/////////////////////////////////////////////////////////////////////// Actual propagation ///////////////////////////////////////////////////////////////////////

/// Perform the integration ///




        Eigen::MatrixXd outputMatrix = performIntegration(desiredOrbitalAltitude,desiredInclinationDeg,initialAltitude,initialLatitudeDeg,initialLongitudeDeg,
                                                          FlightPathAngleDeg,HeadingAngleDeg,initialGroundVelocity,massMAV,thrust,specificImpulse,initialBurnTime,
                                                          constantThrustElevationAngle,constantThrustAzimuthAngle,maxOrder,
                                                          chosenLocalErrorTolerance,chosenStepSize,setEndTime,RKFinitiaterTime,rotatingPlanet,GravityAcc,ThrustAcc,DragAcc,comparison);

    std::cout<<"outputMatrix = "<<outputMatrix<<std::endl;

//////////////////////////////////////////////////////////////////////// Output //////////////////////////////////////////////////////////////////////////////////////////////////
/*    // TSI cartesian end state
    integrationEndStatesAndInfo(0,0) = currentStateAndTime.getCurrentTime();    // end time
    integrationEndStatesAndInfo(0,1) = TSIendState(0);      // x-position
    integrationEndStatesAndInfo(0,2) = TSIendState(1);      // y-position
    integrationEndStatesAndInfo(0,3) = TSIendState(2);      // z-position
    integrationEndStatesAndInfo(0,4) = TSIendState(3);      // x-velocity
    integrationEndStatesAndInfo(0,5) = TSIendState(4);      // y-velocity
    integrationEndStatesAndInfo(0,6) = TSIendState(5);      // z-velocity
    integrationEndStatesAndInfo(0,7) = TSIendState(6);      // MAV mass

    // RKF cartesian end state
    integrationEndStatesAndInfo(1,0) = currentStateAndTime.getCurrentTime();    // end time
    integrationEndStatesAndInfo(1,1) = endState(0);      // x-position
    integrationEndStatesAndInfo(1,2) = endState(1);      // y-position
    integrationEndStatesAndInfo(1,3) = endState(2);      // z-position
    integrationEndStatesAndInfo(1,4) = endState(3);      // x-velocity
    integrationEndStatesAndInfo(1,5) = endState(4);      // y-velocity
    integrationEndStatesAndInfo(1,6) = endState(5);      // z-velocity
    integrationEndStatesAndInfo(1,7) = endState(6);      // MAV mass


    // Final circularisation information
    integrationEndStatesAndInfo(2,0) = deltaVforTSI;    // Delta V required for the circularisation of TSI
    integrationEndStatesAndInfo(2,1) = deltaVforRKF;    // Delta V required for the circularisation of RKF
    integrationEndStatesAndInfo(2,2) = rad2deg(TSIendKeplerElements(2));    // The final inclination of TSI
    integrationEndStatesAndInfo(2,3) = rad2deg(RKFendKeplerElements(2));    // The final inclination of RKF
    integrationEndStatesAndInfo(2,4) = propMassTSI; // The propellant mass used for the circularisation burn of TSI
    integrationEndStatesAndInfo(2,5) = propMassRKF; // The propellant mass used for the circularisation burn of RKF
    integrationEndStatesAndInfo(2,6) = finalMassTSI; // The final total system mass of TSI
    integrationEndStatesAndInfo(2,7) = finalMassRKF; // The final total system mass of RKF

    // Time information
    integrationEndStatesAndInfo(3,0) = tsiTime(0); // TSI CPU time
    integrationEndStatesAndInfo(3,1) = tsiTime(1); // TSI wall time
    integrationEndStatesAndInfo(3,2) = RKFtime(0); // RKF CPU time
    integrationEndStatesAndInfo(3,3) = RKFtime(1); // RKF wall time
    integrationEndStatesAndInfo(3,4) = countTSI; // TSI steps/function evaluations
    integrationEndStatesAndInfo(3,5) = count; // RKF steps evaluations
    integrationEndStatesAndInfo(3,6) = count*13; // RKF function evaluations

    // Difference fraction between TSI and RKF
    integrationEndStatesAndInfo(4,0) = differenceFractionEnd(0); // Difference fraction in x-position
    integrationEndStatesAndInfo(4,1) = differenceFractionEnd(1); // Difference fraction in y-position
    integrationEndStatesAndInfo(4,2) = differenceFractionEnd(2); // Difference fraction in z-position
    integrationEndStatesAndInfo(4,3) = differenceFractionEnd(3); // Difference fraction in x-velocity
    integrationEndStatesAndInfo(4,4) = differenceFractionEnd(4); // Difference fraction in y-velocity
    integrationEndStatesAndInfo(4,5) = differenceFractionEnd(5); // Difference fraction in z-velocity
    integrationEndStatesAndInfo(4,6) = differenceFractionEnd(6); // Difference fraction in mass

    // Difference fraction between TSI and RKF
    integrationEndStatesAndInfo(5,0) = differenceEnd(0); // Difference in x-position
    integrationEndStatesAndInfo(5,1) = differenceEnd(1); // Difference in y-position
    integrationEndStatesAndInfo(5,2) = differenceEnd(2); // Difference in z-position
    integrationEndStatesAndInfo(5,3) = differenceEnd(3); // Difference in x-velocity
    integrationEndStatesAndInfo(5,4) = differenceEnd(4); // Difference in y-velocity
    integrationEndStatesAndInfo(5,5) = differenceEnd(5); // Difference in z-velocity
    integrationEndStatesAndInfo(5,6) = differenceEnd(6); // Difference in mass

//*/    /// Output ///

    /// Storing the output TSI


//    outputMatrixTSI(run,0) = maxOrder;
//    outputMatrixTSI(run,0) = chosenLocalErrorTolerance;
//    outputMatrixTSI(run,0) = initialAltitude;
//        outputMatrixTSI(run,0) = initialLatitudeDeg;
//        outputMatrixTSI(run,0) = FlightPathAngleDeg;
        outputMatrixTSI(run,0) = HeadingAngleDeg;
    //    outputMatrixTSI(run,0) = initialGroundVelocity;
//            outputMatrixTSI(run,0) = run+1;

//        outputMatrixTSI(run,1) = initialLongitudeDeg;

    outputMatrixTSI(run,2) = outputMatrix(3,0); // CPU time TSI
    outputMatrixTSI(run,3) = outputMatrix(3,1); // Wall time TSI
    outputMatrixTSI(run,4) = outputMatrix(3,4); // Number of evaluations TSI

    // State and time
    outputMatrixTSI(run,5) = outputMatrix(0,0); // Time
    outputMatrixTSI(run,6) = outputMatrix(0,1);
    outputMatrixTSI(run,7) = outputMatrix(0,2);
    outputMatrixTSI(run,8) = outputMatrix(0,3);
    outputMatrixTSI(run,9) = outputMatrix(0,4);
    outputMatrixTSI(run,10) = outputMatrix(0,5);
    outputMatrixTSI(run,11) = outputMatrix(0,6);
    outputMatrixTSI(run,12) = outputMatrix(0,7); // Mass

    // Fraction difference w.r.t. RKF
    outputMatrixTSI(run,13) = outputMatrix(4,0);
    outputMatrixTSI(run,14) = outputMatrix(4,1);
    outputMatrixTSI(run,15) = outputMatrix(4,2);
    outputMatrixTSI(run,16) = outputMatrix(4,3);
    outputMatrixTSI(run,17) = outputMatrix(4,4);
    outputMatrixTSI(run,18) = outputMatrix(4,5);
    outputMatrixTSI(run,19) = outputMatrix(4,6);

    // Value difference w.r.t. RKF
    outputMatrixTSI(run,20) = outputMatrix(5,0);
    outputMatrixTSI(run,21) = outputMatrix(5,1);
    outputMatrixTSI(run,22) = outputMatrix(5,2);
    outputMatrixTSI(run,23) = outputMatrix(5,3);
    outputMatrixTSI(run,24) = outputMatrix(5,4);
    outputMatrixTSI(run,25) = outputMatrix(5,5);
    outputMatrixTSI(run,26) = outputMatrix(5,6);

    // In case of circularisation:
    outputMatrixTSI(run,27) = outputMatrix(2,4); // The propellant mass used for the circularisation burn of TSI

    /// Storing the output RKF

    outputMatrixRKF(run,0) = outputMatrixTSI(run,0);
    outputMatrixRKF(run,1) = outputMatrixTSI(run,1);

    outputMatrixRKF(run,2) = outputMatrix(3,2); // CPU time RKF
    outputMatrixRKF(run,3) = outputMatrix(3,3); // Wall time RKF
    outputMatrixRKF(run,4) = outputMatrix(3,6); // Number of function evaluations RKF

    // State and time
    outputMatrixRKF(run,5) = outputMatrix(1,0); // Time
    outputMatrixRKF(run,6) = outputMatrix(1,1);
    outputMatrixRKF(run,7) = outputMatrix(1,2);
    outputMatrixRKF(run,8) = outputMatrix(1,3);
    outputMatrixRKF(run,9) = outputMatrix(1,4);
    outputMatrixRKF(run,10) = outputMatrix(1,5);
    outputMatrixRKF(run,11) = outputMatrix(1,6);
    outputMatrixRKF(run,12) = outputMatrix(1,7); // Mass

    // In case of circularisation:
    outputMatrixRKF(run,13) = outputMatrix(2,5); // The propellant mass used for the circularisation burn of RKF




    run++;


    } // End for
    /////////// End of for loop //////////////



    /// Storing the output in a .csv file
    // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
    const std::string outputDirectoryTSI = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/04.VerificationAndValidation/04.DifferentVariables/" + rotOrNot + "/TSI/" + currentVariable + "/" + caseName + "/";

    // Set output format for matrix output.
    Eigen::IOFormat csvFormatTSI( 15, 0, ", ", "\n" );

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

//        std::cout<<"currentSec = "<<currentSec<<std::endl;

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

    std::string ComputerTimeString = currentYear + "-" + currentMonth + "-" + currentDay + "_" + currentHour + ":" + currentMin;  // Convert to string and store

    std::string newFileName = currentVariable + "DataCollectionFile_TSI_" + caseName + "_" + ComputerTimeString + ".csv";



// Set new absolute path to file containing the data.
std::string dataAbsolutePathTSI = outputDirectoryTSI + newFileName;



    // Storing the data

    std::ifstream ifileTSI(dataAbsolutePathTSI.c_str()); // Check it as an input file

    bool fexistsTSI = false;   // Set the default to "It does not exist"

    if (ifileTSI){         // Attempt to open the file


       fexistsTSI = true;      // If the file can be opened it must exist

       ifileTSI.close();   // Close the file

    }


    // If so: error and create temporary file, if not: create new file and put data in

    if (fexistsTSI == true){

    std::cerr<<"Error, this file name already excists, please use a different file name to store the data, or simply wait one minute before running again"<<std::endl;


}
        else{

        // Export the data.
        std::ofstream exportFile1( dataAbsolutePathTSI.c_str( ) ); // Make the new file
        std::cout<<"New file called "<<dataAbsolutePathTSI<<" has been created"<<std::endl;
        exportFile1 << outputMatrixTSI.format( csvFormatTSI );          // Store the new values
        exportFile1.close( );   // Close the file


    };


    /// Storing the RKF values

    /// Storing the output in a .csv file
    // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
    const std::string outputDirectoryRKF = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/04.VerificationAndValidation/04.DifferentVariables/" + rotOrNot + "/RKF/" + currentVariable + "/" + caseName + "/";

    // Set output format for matrix output.
    Eigen::IOFormat csvFormatRKF( 15, 0, ", ", "\n" );



    std::string newFileNameRKF = currentVariable + "DataCollectionFile_RKF_" + caseName + "_" + ComputerTimeString + ".csv";



// Set new absolute path to file containing the data.
std::string dataAbsolutePathRKF = outputDirectoryRKF + newFileNameRKF;



    // Storing the data

    std::ifstream ifileRKF(dataAbsolutePathRKF.c_str()); // Check it as an input file

    bool fexistsRKF = false;   // Set the default to "It does not exist"

    if (ifileRKF){         // Attempt to open the file


       fexistsRKF = true;      // If the file can be opened it must exist

       ifileRKF.close();   // Close the file

    }


    // If so: error and create temporary file, if not: create new file and put data in

    if (fexistsRKF == true){

    std::cerr<<"Error, this RKF file name already excists, please use a different file name to store the data, or simply wait one minute before running again"<<std::endl;


}
        else{

        // Export the data.
        std::ofstream exportFile1( dataAbsolutePathRKF.c_str( ) ); // Make the new file
        std::cout<<"New file called "<<dataAbsolutePathRKF<<" has been created"<<std::endl;
        exportFile1 << outputMatrixRKF.format( csvFormatRKF );          // Store the new values
        exportFile1.close( );   // Close the file


    };


    return 0;
}


