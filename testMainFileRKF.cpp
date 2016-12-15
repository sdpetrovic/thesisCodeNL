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
 *      160518    S.D. Petrovic     Changed the state derivative function header/source file structure into a class. Now the input for the function itself is simply time and state.
 *                                  This way, it should be able to be used by the existing RKF integrator from the tudat library.
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

// Used by my own functions
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

/// Testing the auxiliary equations ///
#include <thesisProject/Auxiliary.h>                // Original test file

/// Testing the other required functions ///
#include <thesisProject/projectLibraries/otherRequiredFunctions.h>              // Original test file


/* Taylor Series:
 *
/// Testing the basic recurrence relations ///
#include <thesisProject/projectLibraries/basicRecurrenceRelations.h>               // Original test file

/// Testing all recurrence relations ///
#include <thesisProject/projectLibraries/allRecurrenceRelations.h>          // Original test file

/// Testing the stepSize class ///
#include <thesisProject/StepSize.h>             // Original test file

/// Testing the actual Taylor Series integration fucntion ///
#include <thesisProject/projectLibraries/TaylorSeriesIntegration.h>             // Original test file

//*/


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
    const double rotationalVelocityMars = Mars.rotationalVelocity();
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

    const Eigen::Vector3d initialCartesianPositionInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocityMars*inertialFrameTime-primeMeridianAngle)*initialCartesianPositionRotationalFrame;


    // Compute initial velocity in y-direction as seen from the launch site in the inertial frame

    const Eigen::Vector3d initialVelocityLaunchSite = Eigen::Vector3d(0,(rotationalVelocityMars*initialRadius*cos(initialLatitude)),0);

    const Eigen::Vector3d initialVelocityInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocityMars*inertialFrameTime-primeMeridianAngle+initialLongitude)*initialVelocityLaunchSite;

    /// Setting StateAndTime class using modified vector ///

    tudat::basic_mathematics::Vector7d aState;

    aState(0) = initialCartesianPositionInertialFrame(0);       // x,y,z position in [km]
    aState(1) = initialCartesianPositionInertialFrame(1);
    aState(2) = initialCartesianPositionInertialFrame(2);
    aState(3) = initialVelocityInertialFrame(0);                // x,y,z velocity in [km/s]
    aState(4) = initialVelocityInertialFrame(1);
    aState(5) = initialVelocityInertialFrame(2);
    aState(6) = 227;  // Mass [kg] from literature study

    StateAndTime currentStateAndTime(aState);        // Creating the current state class using the namespace and class directly





///////////////////////////////////////////////////////////////////////////////////
////////////////////// Testing the State Derivative Function //////////////////////
///////////////////////////////////////////////////////////////////////////////////

    /// Setting the data collection file for RKF and inserting the first values ///

        // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
        const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/01.integrationResults/RKF/";


        // Set output format for matrix output.
        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

        // Set absolute path to file containing the data.
        std::string dataAbsolutePath = outputDirectory + "test1RKFstateAndTime.csv";

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

            std::string ComputerTimeString = currentHour + ":" + currentMin + ":" + currentSec;  // Convert to string and store

            std::string newFileName = "backupRKFFileAtTime_" + ComputerTimeString + ".csv";

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


    /// Testing with the state derivative function header and source file
//    const tudat::basic_mathematics::Vector7d stateDerivativeVector = ascentStateDerivativeFunction(Mars,MAV,stateAndTime);

//    std::cout<<"The state derivative vector is "<<stateDerivativeVector<<std::endl;


    /// Testing with the state derivative function class

    // Initial test
//    ascentStateDerivativeFunctionClass stateDerivativeFunctionClass;    // Initialize the class

//    const tudat::basic_mathematics::Vector7d stateDerivativeClassVector = stateDerivativeFunctionClass.ascentStateDerivativeFunction(Mars,MAV,stateAndTime);

//    std::cout<<"The state derivative vector obtained from the class is "<<stateDerivativeClassVector<<std::endl;

    // Full complete test

    ascentStateDerivativeFunctionClass stateDerivativeFunctionClass(Mars,MAV);     // Initialize the class

//    const tudat::basic_mathematics::Vector7d stateDerivativeClassVector = stateDerivativeFunctionClass.ascentStateDerivativeFunction(stateAndTime.getCurrentTime(),stateAndTime.getCurrentState());

//    std::cout<<"The state derivative vector obtained from the class is "<<stateDerivativeClassVector<<std::endl;

/// Creating pointer for the state derivative function, such that it can be called from the integrator ///  as by the boost tutorial of the Tudat wiki
/*
    typedef boost::shared_ptr< ascentStateDerivativeFunctionClass > stateDerivativePointer; // Create a typedef for the shared_ptr to an ascentStateDerivativeFunctionClass

    stateDerivativePointer exampleStateDerivativePointer( new ascentStateDerivativeFunctionClass(Mars,MAV)); // Create a boost shared pointer

    stateDerivativePointer newExampleStateDerivativePointer = boost::make_shared< ascentStateDerivativeFunctionClass >(Mars,MAV);   // Create a boost shared pointer using boost::make_shared instead of new
//*/

    // Just using boost::bind

    boost::function< tudat::basic_mathematics::Vector7d( const double, const tudat::basic_mathematics::Vector7d ) > stateDerivativeFunction // Using boost::function to create the passing function
            =boost::bind( &ascentStateDerivativeFunctionClass::ascentStateDerivativeFunction, &stateDerivativeFunctionClass, _1, _2);       // Then using boost::bind to bind the function as per class stateDerivativeFunctionClass which requires two inputs:
                                                                                                                                            // _1 which is the current time and _2 which is the current state
    

 ///// Testing the implementation in the integrator ///

    // Initial conditions.
    const double initialTime = currentStateAndTime.getCurrentTime();                            // Time.
//    Eigen::VectorXd initialState = currentStateAndTime.getCurrentState(); // State: start with zero velocity at the origin.

    const double endTime = 0.2;     // Using the same initial step-size as defined for TSI

    // Step-size settings.
    // The minimum and maximum step-size are set such that the input data is fully accepted by the
    // integrator, to determine the steps to be taken.
    const double zeroMinimumStepSize = std::numeric_limits< double >::epsilon( );
    const double infiniteMaximumStepSize = std::numeric_limits< double >::infinity( );
    double stepSize = 0.2;          // Using the same initial step-size as defined for TSI

    // Tolerances.
    const double relativeTolerance = 1e-3;     // 1e-14 is used by TSI, original setting was 1e-15
    const double absoluteTolerance = 1e-3;     // 1e-14 is used by TSI, original setting was 1e-15

    // For RKF7(8) a step-size of 0.2 is only used if the tolerances are 1e-3.... and is accepted till 1e-8
    // For RKF4(5) a step-size of 0.2 is only used if the tolerances are 1e-7.... and is accepted till 1e-10
    // For DP8(7) a step-size of 0.2 is only used if the tolerances are 1e-7.... and is accepted till 1e-9

    /// RungeKutta4 numerical integrator.

    tudat::numerical_integrators::RungeKutta4IntegratorXd RK4integrator(
        stateDerivativeFunction, initialTime, initialState );

    // Integrate to the specified end time.
    Eigen::VectorXd RK4endState = RK4integrator.integrateTo( endTime, stepSize );

    std::cout<<"The RK4 end state is "<<RK4endState<<std::endl;

    /// Runge-Kutta-Fehlberg 7(8) integrator.
       tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
                   tudat::numerical_integrators::RungeKuttaCoefficients::get(
                       tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78),
                   stateDerivativeFunction, initialTime, initialState, zeroMinimumStepSize,
                   infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

//       /// Runge-Kutta-Fehlberg 4(5) integrator.
//          tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
//                      tudat::numerical_integrators::RungeKuttaCoefficients::get(
//                          tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg45),
//                      stateDerivativeFunction, initialTime, initialState, zeroMinimumStepSize,
//                      infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

//          /// Dormand-Prince 8(7) integrator.
//             tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
//                         tudat::numerical_integrators::RungeKuttaCoefficients::get(
//                             tudat::numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince),
//                         stateDerivativeFunction, initialTime, initialState, zeroMinimumStepSize,
//                         infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

/*/// Debug ///
//          std::cout<<"They initialize just fine..."<<std::endl;

       // Perform a single integration step.
       integrator.performIntegrationStep( stepSize );



       // The result of the integration.
       Eigen::VectorXd RKFendState = integrator.getCurrentState( );

       std::cout<<"The RKF end state is "<<RKFendState<<std::endl;
       std::cout<<"The running time = "<<integrator.getCurrentIndependentVariable()<<std::endl;
       std::cout<<"The next step-size = "<<integrator.getNextStepSize()<<std::endl;
/// Debug ///
       //*/

       // Set initial running time. This is updated after each step that the numerical integrator takes.
    double runningTime = 0.0;
    int count = 0;

    do
    {
        // Make sure the integrator does not integrate beyond the end time.



        if ( std::fabs( endTime - runningTime )
             <= std::fabs( stepSize ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) )
        {
            stepSize = endTime - integrator.getCurrentIndependentVariable( );
        }

        double prevStepSize = stepSize;

        // Perform a single integration step. Then update the step-size and running time.
        integrator.performIntegrationStep( stepSize );
        stepSize = integrator.getNextStepSize( );
        runningTime = integrator.getCurrentIndependentVariable( );

        count++;

        Eigen::VectorXd currentState = integrator.getCurrentState();

        std::cout<<"The current stepSize is "<<prevStepSize<<" s"<<std::endl;


        /// Storing the values ///

                outputVector = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.

                // Filling the output vector
                outputVector(0,0) = integrator.getCurrentIndependentVariable();   // Storing the updated time
                outputVector(0,1) = currentState(0);   // Storing the updated x position
                outputVector(0,2) = currentState(1);   // Storing the updated y position
                outputVector(0,3) = currentState(2);   // Storing the updated z position
                outputVector(0,4) = currentState(3);   // Storing the updated x velocity
                outputVector(0,5) = currentState(4);   // Storing the updated y velocity
                outputVector(0,6) = currentState(5);   // Storing the updated z velocity
                outputVector(0,7) = currentState(6);   // Storing the updated MAV mass


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






    }while( !( endTime - runningTime <= std::numeric_limits< double >::epsilon( ) ) );

    // The result of the integration.
    Eigen::VectorXd endState = integrator.getCurrentState( );

    std::cout<<"Final number of integration steps is "<<count<<std::endl;
    std::cout<<"The end state is "<<endState<<std::endl;




    return 0;
}


