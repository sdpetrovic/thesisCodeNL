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
 *      160411    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

// This is a test main file to test the different class files, header/source files to see if any output is produced (verification)

#include <iostream>
#include <Eigen/Core>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <string>

#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

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

#include <thesisProject/stateAndTime.h>     // Original test file


int main()

{

/*

    /// Testing the Celestial Body class ///

//    // First test

//    celestialBody Mars;

    // Second test

//    const std::string planet = "Mars";
//    const std::string planet = "Venus";

    celestialBody Mars;

//    celestialBody Mars(planet);
//    Mars.setPlanet(planet);

   const double adiabeticIndex = Mars.adiabeticIndex();
   const double specificGasConstant = Mars.specificGasConstant();
   const double standardGravitationalParameter = Mars.standardGravitationalParameter();
   const double rotationalVelocity = Mars.rotationalVelocity();
   const double primeMeridianAngle = Mars.primeMeridianAngle();
   const double inertialFrameTime = Mars.inertialFrameTime();

   const Eigen::MatrixXd temperaturePolyCoefficients = Mars.temperaturePolyCoefficients();
   const Eigen::MatrixXd temperatureAltitudeRanges = Mars.temperatureAltitudeRanges();
   const Eigen::VectorXd densityPolyCoefficients = Mars.densityPolyCoefficients();

       std::cout<<"The adiabeticIndex for the Martian atmosphere is "<<adiabeticIndex<<std::endl;
       std::cout<<"The specificGasConstant for the Martian atmosphere is "<<specificGasConstant<<std::endl;
       std::cout<<"The standardGravitationalParameter for the Martian atmosphere is "<<standardGravitationalParameter<<std::endl;
       std::cout<<"The rotationalVelocity for the Martian atmosphere is "<<rotationalVelocity<<std::endl;
       std::cout<<"The primeMeridianAngle for the Martian atmosphere is "<<primeMeridianAngle<<std::endl;
       std::cout<<"The inertialFrameTime for the Martian atmosphere is "<<inertialFrameTime<<std::endl;
       std::cout<<"The temperaturePolyCoefficients for the Martian atmosphere is "<<temperaturePolyCoefficients<<std::endl;
       std::cout<<"The temperatureAltitudeRanges for the Martian atmosphere is "<<temperatureAltitudeRanges<<std::endl;
       std::cout<<"The densityPolyCoefficients for the Martian atmosphere is "<<densityPolyCoefficients<<std::endl;



//   boost::shared_ptr< celestialBody > MarsTwo = boost::make_shared< celestialBody > ();     // Still not sure how this works so lets keep it simple for now shall we...

//     boost::shared_ptr< Mars.adiabeticIndex() > adiabeticIndexTwoTest = boost::make_shared< Mars.adiabeticIndex() > ();

//            std::cout<<"The adiabeticIndexTwoTest for the Martian atmosphere is "<<adiabeticIndexTwoTest<<std::endl;


    /// Testing the vehicle class ///

    MarsAscentVehicle MAV;

    // Returning the different constant parameters

    const double Thrust = MAV.Thrust();                                                          // T     engine nominal thrust
    const double specificImpulse = MAV.specificImpulse();                                        // Isp     engine nominal specific impulse
    const double referenceArea = MAV.referenceArea();                                             // S   vehicle reference area
    const Eigen::MatrixXd dragCoefficientPolyCoefficients = MAV.dragCoefficientPolyCoefficients();            // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
    const Eigen::MatrixXd dragCoefficientMachranges = MAV.dragCoefficientMachranges();                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient
    const Eigen::MatrixXd thrustAzimuth = MAV.thrustAzimuth();                                // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
    const Eigen::MatrixXd thrustElevation = MAV.thrustElevation();                                // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)

    std::cout<<"The Thrust for the Mars Ascent Vehicle is "<<Thrust<<std::endl;
    std::cout<<"The specificImpulse for the Mars Ascent Vehicle is "<<specificImpulse<<std::endl;
    std::cout<<"The referenceArea for the Mars Ascent Vehicle is "<<referenceArea<<std::endl;
    std::cout<<"The dragCoefficientPolyCoefficients for the Mars Ascent Vehicle is "<<dragCoefficientPolyCoefficients<<std::endl;
    std::cout<<"The dragCoefficientMachranges for the Mars Ascent Vehicle is "<<dragCoefficientMachranges<<std::endl;
    std::cout<<"The thrustAzimuth for the Mars Ascent Vehicle is "<<thrustAzimuth<<std::endl;
    std::cout<<"The thrustElevation for the Mars Ascent Vehicle is "<<thrustElevation<<std::endl;


    /// Testing the set functions for the MAV ///

    // Thrust Azimuth-Gimbal Angles
       Eigen::MatrixXd thrustAzimuth_ = Eigen::MatrixXd::Zero(6,3); // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)

        // Section 1
        thrustAzimuth_(0,0) = 0.0;   // Lower bound time
        thrustAzimuth_(0,1) = 1.0;   // Upper bound time

        thrustAzimuth_(0,2) = 0.0;   // Thrust azimuth angle

        // Section 2
        thrustAzimuth_(1,0) = 1.0;   // Lower bound time
        thrustAzimuth_(1,1) = 10.0;   // Upper bound time

        thrustAzimuth_(1,2) = 0.05;   // Thrust azimuth angle

        // Section 3
        thrustAzimuth_(2,0) = 10.0;   // Lower bound time
        thrustAzimuth_(2,1) = 100.0;   // Upper bound time

        thrustAzimuth_(2,2) = 0.176;   // Thrust azimuth angle

        // Section 4
        thrustAzimuth_(3,0) = 100.0;   // Lower bound time
        thrustAzimuth_(3,1) = 120.0;   // Upper bound time

        thrustAzimuth_(3,2) = 0.1;   // Thrust azimuth angle

        // Section 5
        thrustAzimuth_(4,0) = 120.0;   // Lower bound time
        thrustAzimuth_(4,1) = 200.0;   // Upper bound time

        thrustAzimuth_(4,2) = 0.5;   // Thrust azimuth angle

        // Section 6
        thrustAzimuth_(5,0) = 200.0;   // Lower bound time
        thrustAzimuth_(5,1) = 400.0;   // Upper bound time

        thrustAzimuth_(5,2) = 0.0;   // Thrust azimuth angle


        // Thrust Elevation-Gimbal Angles
         Eigen::MatrixXd   thrustElevation_ = Eigen::MatrixXd::Zero(6,3); // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)

            // Section 1
            thrustElevation_(0,0) = 0.0;   // Lower bound time
            thrustElevation_(0,1) = 1.0;   // Upper bound time

            thrustElevation_(0,2) = 0.0;   // Thrust elevation angle

            // Section 2
            thrustElevation_(1,0) = 1.0;   // Lower bound time
            thrustElevation_(1,1) = 10.0;   // Upper bound time

            thrustElevation_(1,2) = 0.05;   // Thrust elevation angle

            // Section 3
            thrustElevation_(2,0) = 10.0;   // Lower bound time
            thrustElevation_(2,1) = 100.0;   // Upper bound time

            thrustElevation_(2,2) = 0.176;   // Thrust elevation angle

            // Section 4
            thrustElevation_(3,0) = 100.0;   // Lower bound time
            thrustElevation_(3,1) = 120.0;   // Upper bound time

            thrustElevation_(3,2) = 0.1;   // Thrust elevation angle

            // Section 5
            thrustElevation_(4,0) = 120.0;   // Lower bound time
            thrustElevation_(4,1) = 200.0;   // Upper bound time

            thrustElevation_(4,2) = 0.5;   // Thrust elevation angle

            // Section 6
            thrustElevation_(5,0) = 200.0;   // Lower bound time
            thrustElevation_(5,1) = 400.0;   // Upper bound time

            thrustElevation_(5,2) = 0.0;   // Thrust elevation angle





      // Set the new matrices and test it

        MAV.setThrustAzimuth(thrustAzimuth_);

        const Eigen::MatrixXd thrustAzimuthTwo = MAV.thrustAzimuth();                                // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
        const Eigen::MatrixXd thrustElevationTwo = MAV.thrustElevation();                                // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)

        std::cout<<"The thrustAzimuth for the Mars Ascent Vehicle is "<<thrustAzimuthTwo<<std::endl;
        std::cout<<"The thrustElevation for the Mars Ascent Vehicle is "<<thrustElevationTwo<<std::endl;

        MAV.setThrustElevation(thrustElevation_);

        const Eigen::MatrixXd thrustAzimuthThree = MAV.thrustAzimuth();                                // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
        const Eigen::MatrixXd thrustElevationThree = MAV.thrustElevation();                                // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)

        std::cout<<"The thrustAzimuth for the Mars Ascent Vehicle is "<<thrustAzimuthThree<<std::endl;
        std::cout<<"The thrustElevation for the Mars Ascent Vehicle is "<<thrustElevationThree<<std::endl;
*/

/// Testing the current state and time class ///


    // Initial conditions

    tudat::basic_mathematics::Vector6d aState;

    aState(0,0) = 1;
    aState(0,1) = 1;
    aState(0,2) = 1;
    aState(0,3) = 2;
    aState(0,4) = 2;
    aState(0,5) = 2;

    double aMass = 227;   // [kg] from literature study

    satellite_propagator_examples::StateAndTime currentStateAndTime(aState,aMass);




    return 0;
}
