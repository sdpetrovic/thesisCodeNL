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
 *      160505    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

#include "TaylorSeriesIntegration.h"



/* ///// Some required functions ///// (Moved to separate header and source file)

/// deg2rad and rad2deg ///

const double deg2rad(const double deg){

    const double rad = deg*tudat::mathematical_constants::LONG_PI/180;

    return rad;
}

const double rad2deg(const double rad){

    const double deg = rad*180/tudat::mathematical_constants::LONG_PI;

    return deg;
}


/// B-P frame transformations ///



//! Get transformation quaternion from the Body (B) to the Propulsion (P) frame.
Eigen::Quaterniond getBodyToPropulsionFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd RotationAroundZaxis = Eigen::AngleAxisd(
                -1.0 * thrustAzimuth, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd RotationAroundYaxis = Eigen::AngleAxisd(
                -1.0 * thrustElevation,
                Eigen::Vector3d::UnitY( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( RotationAroundYaxis * RotationAroundZaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}


//! Get transformation matrix from the Body (B) to the Propulsion (P) frame.
Eigen::Matrix3d getBodyToPropulsionFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation )
{
    return getBodyToPropulsionFrameTransformationQuaternion(
            thrustAzimuth, thrustElevation ).toRotationMatrix( );
}

//! Get transformation matrix from the Propulsion (P) to the Body (B) frame.
Eigen::Matrix3d getPropulsionToBodyFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation )
{
    return getBodyToPropulsionFrameTransformationMatrix(
            thrustAzimuth, thrustElevation ).transpose( );
}

//! Get transformation quaternion from the Propulsion (P) to the Body (B) frame.
Eigen::Quaterniond getPropulsionToBodyFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation )
{
    return getBodyToPropulsionFrameTransformationQuaternion(
            thrustAzimuth, thrustElevation ).inverse( );
}
//*/

/// Taylor series integration step ///


Eigen::VectorXd performTaylorSeriesIntegrationStep(const celestialBody& planet_, const MarsAscentVehicle& MAV_, const StateAndTime& currentStateAndTime_, StepSize& stepSize, const double maxOrder_ ,
                                                   const double FlightPathAngle_,
                                                   const double HeadingAngle_){

    // The stepSize class is the only class that will be updated directly from this function!

        // Create variables to be used in this function

    celestialBody planet = planet_;                               // The celestial body class
    MarsAscentVehicle MAV = MAV_;                                 // The MAV class
    StateAndTime currentStateAndTime = currentStateAndTime_;      // The current state and time class
    const double maxOrder = maxOrder_;                                  // The maximum order K (default value is 20)

    // Getting the required information from the classes

    // celestialBody
    const double adiabeticIndex = planet.adiabeticIndex();                  // The adiabetic index gamma_a
    const double specificGasConstant = planet.specificGasConstant();        // The specific gas constant Rstar
    const double standardGravitationalParameter = planet.standardGravitationalParameter();  // The standard gravitational parameter mu
    const double rotationalVelocity = planet.rotationalVelocity();          // The rotational velocity of the planet Omega
    const double primeMeridianAngle = planet.primeMeridianAngle();          // The angle of the prime meridian when the inertial frame was set omega_p
    const double inertialFrameTime = planet.inertialFrameTime();            // Time between start of simulation and when inertial frame was set

    const double bodyReferenceRadius = planet.bodyReferenceRadius();         // Reference radius of the planet R


    const Eigen::MatrixXd temperaturePolyCoefficients = planet.temperaturePolyCoefficients();   // Temperature polynomial coefficients for temperature curve P_Tn
    const Eigen::MatrixXd temperatureAltitudeRanges = planet.temperatureAltitudeRanges();       // Temperature curve altitude ranges h
    const Eigen::VectorXd densityPolyCoefficients = planet.densityPolyCoefficients();           // Density polynomial coefficients for density curve P_rho n

    // MarsAscentVehicle
    const double Thrust = MAV.Thrust();                                                          // T     engine nominal thrust
    const double specificImpulse = MAV.specificImpulse();                                        // Isp     engine nominal specific impulse
    const double referenceArea = MAV.referenceArea();                                             // S   vehicle reference area
    const Eigen::MatrixXd dragCoefficientPolyCoefficients = MAV.dragCoefficientPolyCoefficients();            // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
    const Eigen::MatrixXd dragCoefficientMachRanges = MAV.dragCoefficientMachRanges();                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient
    const Eigen::MatrixXd thrustAzimuth = MAV.thrustAzimuth();                                // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
    const Eigen::MatrixXd thrustElevation = MAV.thrustElevation();                                // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)

    // StateAndTime

    const tudat::basic_mathematics::Vector7d currentState = currentStateAndTime.getCurrentState(); // The complete state including position, velocity and mass
    const double currentMass = currentStateAndTime.getCurrentMass();                                // The mass seperately
    const double currentTime = currentStateAndTime.getCurrentTime();                                // The current time

    const double currentStepSize = stepSize.getCurrentStepSize();                      // The current step-size

    // Specified initial conditions

    const double FlightPathAngle = FlightPathAngle_;            // Set flight-path angle
    const double HeadingAngle = HeadingAngle_;                  // Set heading angle

//    std::cout<<"Flight-path angle in integration.cpp = "<<FlightPathAngle<<std::endl;
//    std::cout<<"Heading angle in integration.cpp = "<<HeadingAngle<<std::endl;


    //std::cout<<"This works right 1?"<<std::endl;

//////// Computations //////////

    /// Thrust acceleration in B-frame ///   thrustAccelerationsBframe // initial (for u)

    const Eigen::Vector3d thrustAccelerationsPframe = Eigen::Vector3d((Thrust/currentMass),0,0);            // THIS HAS TO BE CHANGED IN THE FUTURE TO INCLUDE A WIDE RANGE OF THRUST AZIMUTH AND ELEVATION ANGLES!!!

    Eigen::MatrixXd thrustAzimuthMatrix = MAV.thrustAzimuth();
    Eigen::MatrixXd thrustElevationMatrix = MAV.thrustElevation();

//    const double thrustAzimuthTestDeg = 0;             // thrust azimuth gimbal angle [Deg] 10 for testing
//    const double thrustElevationTestDeg = 0;            // thrust elevation gimbal angle [Deg] 5 for testing

//    const double thrustAzimuthTest = deg2rad(thrustAzimuthTestDeg);     // thrust azimuth gimbal angle [rad]
//    const double thrustElevationTest = deg2rad(thrustElevationTestDeg); // thrust elevation gimbal angle [rad]

      const double thrustAzimuthTest = thrustAzimuthMatrix(0,2);             // thrust azimuth gimbal angle [rad]
      const double thrustElevationTest = thrustElevationMatrix(0,2);            // thrust elevation gimbal angle [rad]


    const Eigen::Vector3d thrustAccelerationsBframe = getPropulsionToBodyFrameTransformationMatrix(thrustAzimuthTest,thrustElevationTest)*thrustAccelerationsPframe;

//    std::cout<<"thrustAccelerationsBframe = "<<thrustAccelerationsBframe<<std::endl;

//std::cout<<"This works right 2?"<<std::endl;

    /// Computing the auxiliary equations, derivatives and functions ///

    Auxiliary Aux(adiabeticIndex, specificGasConstant,standardGravitationalParameter, rotationalVelocity, primeMeridianAngle,
              inertialFrameTime, bodyReferenceRadius, temperaturePolyCoefficients, temperatureAltitudeRanges,
              densityPolyCoefficients, Thrust, thrustAzimuthMatrix, thrustElevationMatrix, specificImpulse,
              referenceArea, dragCoefficientPolyCoefficients, dragCoefficientMachRanges);

    //std::cout<<"This works right 3?"<<std::endl;

    // Set the initial launch angles
    Aux.setFlightPathAngleAndHeadingAngle(FlightPathAngle,HeadingAngle);    // Used to specify a certain flight-path angle other than 90 degrees and heading angle other than 0 degrees at t = 0 sec

    // Compute the auxiliary equations

    Eigen::VectorXd auxiliaryEquations =  Aux.getAuxiliaryEquations(currentState,currentTime,thrustAccelerationsBframe);

//    std::cout<<"The auxiliaryEquations are "<<auxiliaryEquations<<std::endl;

    /// Debug ///
//    if (currentTime == 0.2){
//        auxiliaryEquations(23) = 1;
//    }
    /// Debug ///

//std::cout<<"This works right 4?"<<std::endl;

    // Compute the auxiliary derivatives

    Eigen::VectorXd auxiliaryDerivatives = Aux.getAuxiliaryDerivatives(currentState,currentTime,thrustAccelerationsBframe,auxiliaryEquations);

//    std::cout<<"The auxiliaryDerivatives are "<<auxiliaryDerivatives<<std::endl;

    // Compute the auxiliary functions

    Eigen::MatrixXd auxiliaryFunctions = Aux.getAuxiliaryFunctions(currentState,currentTime,thrustAccelerationsBframe,auxiliaryEquations,auxiliaryDerivatives);

//    std::cout<<"The auxiliaryFunctions are "<<auxiliaryFunctions<<std::endl;

    /// Computing the Taylor Coefficients ///

    Eigen::MatrixXd TaylorCoefficients = getTaylorCoefficients(adiabeticIndex, specificGasConstant, standardGravitationalParameter, rotationalVelocity, primeMeridianAngle,
                          inertialFrameTime, bodyReferenceRadius,temperaturePolyCoefficients, temperatureAltitudeRanges,
                          densityPolyCoefficients, Thrust, specificImpulse,
                          referenceArea, dragCoefficientPolyCoefficients, dragCoefficientMachRanges,
            thrustAccelerationsBframe,
            thrustAzimuthMatrix,
            thrustElevationMatrix,
            auxiliaryEquations,
            auxiliaryDerivatives,
            auxiliaryFunctions,
            currentTime,
            maxOrder);

    /// Storing the Taylor Coefficients to a file ///

    //std::cout<<"Does this even work1?"<<std::endl;

        Eigen::MatrixXd TaylorCoefficientsOutputMatrix = Eigen::MatrixXd::Zero(7,maxOrder+1);       // Create an output matrix for the file without the first empty row

        TaylorCoefficientsOutputMatrix.row(0) = TaylorCoefficients.row(1);                  // The first line entries are the maxOrder+1 Taylor Series Coefficients for     the position in the x-direction
        TaylorCoefficientsOutputMatrix.row(1) = TaylorCoefficients.row(2);                  // The second line entries are the maxOrder+1 Taylor Series Coefficients for    the position in the y-direction
        TaylorCoefficientsOutputMatrix.row(2) = TaylorCoefficients.row(3);                  // The third line entries are the maxOrder+1 Taylor Series Coefficients for     the position in the z-direction
        TaylorCoefficientsOutputMatrix.row(3) = TaylorCoefficients.row(4);                  // The fourth line entries are the maxOrder+1 Taylor Series Coefficients for    the velocity in the x-direction
        TaylorCoefficientsOutputMatrix.row(4) = TaylorCoefficients.row(5);                  // The fifth line entries are the maxOrder+1 Taylor Series Coefficients for     the velocity in the y-direction
        TaylorCoefficientsOutputMatrix.row(5) = TaylorCoefficients.row(6);                  // The sixth line entries are the maxOrder+1 Taylor Series Coefficients for     the velocity in the z-direction
        TaylorCoefficientsOutputMatrix.row(6) = TaylorCoefficients.row(7);                  // The seventh line entries are the maxOrder+1 Taylor Series Coefficients for   the mass


        /// Start debug ///

        //std::cout<<"Does this even work2? :S"<<std::endl;

/*        std::cout<<"x-position coefficients are: "<<TaylorCoefficientsOutputMatrix.row(0)<<std::endl;
        std::cout<<"y-position coefficients are: "<<TaylorCoefficientsOutputMatrix.row(1)<<std::endl;
        std::cout<<"z-position coefficients are: "<<TaylorCoefficientsOutputMatrix.row(2)<<std::endl;
        std::cout<<"x-velocity coefficients are: "<<TaylorCoefficientsOutputMatrix.row(3)<<std::endl;
        std::cout<<"y-velocity coefficients are: "<<TaylorCoefficientsOutputMatrix.row(4)<<std::endl;
        std::cout<<"z-velocity coefficients are: "<<TaylorCoefficientsOutputMatrix.row(5)<<std::endl;
        std::cout<<"mass coefficients are: "<<TaylorCoefficientsOutputMatrix.row(6)<<std::endl;

 //*/       /// End debug ///

        // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
        const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/02.taylorSeriesCoefficientsOutputFolder/";


        // Set output format for matrix output.
        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

        // Set absolute path to file containing the Taylor Series Coefficients.
        const std::string taylorSeriesCoefficientsAbsolutePath = outputDirectory + "test1TaylorSeriesCoefficients(bugSearch31-05-2016).csv";


        // Check if the file already exists.


        std::ifstream ifile(taylorSeriesCoefficientsAbsolutePath.c_str()); // Check it as an input file

        bool fexists = false;   // Set the default to "It does not exist"

        if (ifile){         // Attempt to open the file


           fexists = true;      // If the file can be opened it must exist

           ifile.close();   // Close the file

        }


        // If so: append, if not: create new file and put data in

        if (fexists == true){

            // Export the Taylor Series Coefficients matrix.
            std::ofstream exportFile1;                          // Define the file as an output file


            exportFile1.open(taylorSeriesCoefficientsAbsolutePath.c_str(),std::ios_base::app);      // Open the file in append mode

            exportFile1 << "\n";                                            // Make sure the new matrix start on a new line

            exportFile1 << TaylorCoefficientsOutputMatrix.format( csvFormat ); // Add the new values

            std::cout<<"The file called "<<taylorSeriesCoefficientsAbsolutePath<<" has been appended"<<std::endl;


            exportFile1.close( );   // Close the file
}
            else{

            // Export the Taylor Series Coefficients matrix.
            std::ofstream exportFile1( taylorSeriesCoefficientsAbsolutePath.c_str( ) ); // Make the new file
            std::cout<<"New file called "<<taylorSeriesCoefficientsAbsolutePath<<" has been created"<<std::endl;
            exportFile1 << TaylorCoefficientsOutputMatrix.format( csvFormat );          // Store the new values
            exportFile1.close( );   // Close the file
        };

//std::cout<<"Does this even work3?"<<std::endl;

        /// Performing the actual Taylor Series expansion for every state variable ///



        tudat::basic_mathematics::Vector7d updatedState = tudat::basic_mathematics::Vector7d::Zero();        // Create a vector for the updatedState and setting it to zero

        for (int n = 0; n<updatedState.size();n++){                 // All variables

        for (int k = 0; k<maxOrder+1;k++){                      // Taylor series summation

            updatedState(n) += TaylorCoefficients((n+1),k)*pow(currentStepSize,k);      // Perform one step of the taylor series expansion and then add it to the previous step

//            if (n==6){
//           std::cout<<"updatedState(6) = "<<updatedState(6)<<std::endl;
//           std::cout<<"currentStepSize = "<<currentStepSize<<std::endl;
//           std::cout<<"TaylorCoefficients((n+1),k) = "<<TaylorCoefficients((n+1),k)<<std::endl;

//            }

        } // Taylor series summation

}   // All variables


//std::cout<<"Does this even work4?"<<std::endl;
        double updatedTime = currentTime+currentStepSize;           // Create the updated time variable


        /// Updating the step-size ///

        tudat::basic_mathematics::Vector7d penultimateCoefficients;         // Define the Xn(K-1) vector
        tudat::basic_mathematics::Vector7d lastCoefficients;                // Define the Xn(K) vector

        for (int i = 0; i<7; i++){                                      // Fill the vectors

            penultimateCoefficients(i)= TaylorCoefficients((i+1),maxOrder-1);       // Xn(K-1)
            lastCoefficients(i)= TaylorCoefficients((i+1),(maxOrder));              // Xn(K)
        }

//        std::cout<<"TaylorSeriesIntegration works till here 1"<<std::endl;

        stepSize.updateStepSizeUsingIteration(penultimateCoefficients, lastCoefficients, maxOrder); // Determining the new step-size and updating the current step-size to that

//        std::cout<<"TaylorSeriesIntegration works till here 2"<<std::endl;

        // Please note that instead of updateStepSizeUsingInteration, you can also use updateStepSizeUsingPreviousStepSize, however, this method does not work if the maximum truncation error estimate is larger
        // than the local error tolerance and requires the previous step-size

        //std::cout<<"Does this even work5?"<<std::endl;

        /// Setting the ouput vector ///

        Eigen::VectorXd updatedStateAndTime(8);

        updatedStateAndTime(0) = updatedState(0);   // Updated position in the x-direction
        updatedStateAndTime(1) = updatedState(1);   // Updated position in the y-direction
        updatedStateAndTime(2) = updatedState(2);   // Updated position in the z-direction
        updatedStateAndTime(3) = updatedState(3);   // Updated velocity in the x-direction
        updatedStateAndTime(4) = updatedState(4);   // Updated velocity in the y-direction
        updatedStateAndTime(5) = updatedState(5);   // Updated velocity in the z-direction
        updatedStateAndTime(6) = updatedState(6);   // Updated mass
        updatedStateAndTime(7) = updatedTime;   // Updated time

//        std::cout<<"TaylorSeriesIntegration works till here 3"<<std::endl;

    return updatedStateAndTime;

} // End of the function
