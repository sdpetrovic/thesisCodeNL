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
 *      160719    S.D. Petrovic     File created
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




int main()

{



std::cout<<setprecision(15)<<"Setting output precision to 15"<<std::endl;

    /// Setting the desired end orbit ///

    double desiredOrbitalAltitude = 320.0; // Desired orbital altitude in kmMOLA (320 km is default)
//    const double desiredEccentricity = 0.0; // Desired orbital eccentricity
    double desiredInclination = deg2rad(45.0); // Desired orbital inclination (45 is default)

    /// Setting the Celestial Body class ///


//    // First test

//    celestialBody Mars;

    // Second test

//    const std::string planet = "Mars";
//    const std::string planet = "Venus";

    celestialBody Mars;

//    Mars.setRotationalVelocity(0); // Set Mars as a non-rotating planet for verification


//    const double adiabeticIndex = Mars.adiabeticIndex();
//    const double specificGasConstant = Mars.specificGasConstant();
//    const double standardGravitationalParameter = Mars.standardGravitationalParameter();

    const double rotationalVelocityMars = Mars.rotationalVelocity();
    std::cout<<"rotationalVelocityMars = "<<rotationalVelocityMars<<std::endl;
    const double primeMeridianAngle = Mars.primeMeridianAngle();
    const double inertialFrameTime = Mars.inertialFrameTime();

   const double bodyReferenceRadius = Mars.bodyReferenceRadius();

//const __float128 pi = 4.0*atan(1);

//std::cout<<"pi = "<<pi<<std::endl;


    /// Setting the vehicle class ///

    MarsAscentVehicle MAV;

    /// Selecting the different accelerations ///

    // No Gravity

//    Mars.setStandardGravitationalParameter(0);
//    Mars.setStandardGravitationalParameter(-7.088e-5);
//    std::cout<<"mu_M = "<<Mars.standardGravitationalParameter()<<std::endl;

    if (Mars.standardGravitationalParameter() == 0){
        std::cout<<"NO GRAVITY"<<std::endl;
    }

    // No Drag

//    MAV.setReferenceArea(0);

    if (MAV.referenceArea() == 0){
        std::cout<<"NO DRAG"<<std::endl;
    }

    // No Thrust

//    MAV.setThrust(0);

    if (MAV.Thrust() == 0){
        std::cout<<"NO THRUST"<<std::endl;
    }


    //////////////////////////////////////////////////////// Test Cases ////////////////////////////////////////////////////////

//    Test case from Woolley 2015 (case 10 SSTO)
    std::cout<<"Test case 1: Woolley 2015 SSTO"<<std::endl;
    desiredOrbitalAltitude = 390.0; // Desired orbital altitude in kmMOLA (320 km is default)
    desiredInclination = deg2rad(45.0); // Desired orbital inclination (45 is default)

    MAV.setMAVmass(267.4);                      // Set the MAV GLOM in [kg]
    MAV.setThrust(3.56);                        // Set the MAV thrust in [kN]
    MAV.setThrustResetValue(MAV.Thrust());      // Set the reset value equal to the original given thrust
    MAV.setSpecificImpulse(256);                // Set the MAV specific impulse [s]
    const double initialBurnTime = 142.5;       // Set the burn time from launch till coast [s]
//    const double burnOutAngle = deg2rad(6.0);  // Set the burn out angle (flight-path angle at end of first burn) [rad]
//    const double finalBurnOutMass = 60.7;       // Set the final burn out mass (empty mass + OS mass + excess propellant mass) [kg]
    const double initialLongitudeDeg = 74.5;     // Set the launch latitude in [deg] (tau)
    const double initialLatitudeDeg = 0.0;        // Starting latitude [deg] initial condition in (delta)
    const double HeadingAngle = deg2rad(90.0);  // Set the launch azimuth [rad] (psi)
    Mars.setUpperAltitudeBound(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the temperature
    MAV.setUpdatedFinalAltitude(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the thrust angles
    const double initialAltitude = -0.6;          // Starting altitude [km MOLA] initial condition is -0.6 km MOLA
    const double initialGroundVelocity = 0.00001;          // Starting velocity in km/s (is suppose to be 0.0...) 0.00001 default at initial step-size of 0.01 sec




//// Test case from Joel (hybrid) case7_3_2016_v33
// std::cout<<"Test case 2: Joel (hybrid) case7_3_2016_v33"<<std::endl;
//    desiredOrbitalAltitude = 3875.19000000064-bodyReferenceRadius; // Desired orbital altitude in kmMOLA (320 km is default)
//    desiredInclination = deg2rad(92.6899999999988); // Desired orbital inclination (45 is default)

//    MAV.setMAVmass(288.95985303149);                      // Set the MAV GLOM in [kg]
//    MAV.setThrust(6.01868886452604);                        // Set the MAV thrust in [kN]
//    MAV.setThrustResetValue(MAV.Thrust());      // Set the reset value equal to the original given thrust
//    MAV.setSpecificImpulse(315.9);                // Set the MAV specific impulse [s]
//    const double initialBurnTime = 99.361911852794;       // Set the burn time from launch till coast [s]
////    const double burnOutAngle = deg2rad(6.0);  // Set the burn out angle (flight-path angle at end of first burn) [rad]
////    const double finalBurnOutMass = 60.7;       // Set the final burn out mass (empty mass + OS mass + excess propellant mass) [kg]
//    const double initialLongitudeDeg = 90.0;     // Set the launch latitude in [deg] (tau)
//    const double initialLatitudeDeg = 0.0;        // Starting latitude [deg] initial condition in (delta)
//    const double HeadingAngle = deg2rad(90.0);  // Set the launch azimuth [rad] (psi)
//    Mars.setUpperAltitudeBound(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the temperature
//    MAV.setUpdatedFinalAltitude(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the thrust angles
//    const double initialAltitude = -0.6;          // Starting altitude [km MOLA] initial condition is -0.6 km MOLA
//    const double initialGroundVelocity = 0.000001;          // Starting velocity in km/s (is suppose to be 0.0...) 0.00001 default at initial step-size of 0.01 sec


//// Test case old verification
//    desiredOrbitalAltitude = 320.0; // Desired orbital altitude in kmMOLA (320 km is default)
//    desiredInclination = deg2rad(45.0); // Desired orbital inclination (45 is default)

////    MAV.setMAVmass(267.4);                      // Set the MAV GLOM in [kg]
////    MAV.setThrust(3.56);                        // Set the MAV thrust in [kN]
////    MAV.setThrustResetValue(MAV.Thrust());      // Set the reset value equal to the original given thrust
////    MAV.setSpecificImpulse(256);                // Set the MAV specific impulse [s]
//    const double initialBurnTime = 142.5;       // Set the burn time from launch till coast [s]
////    const double burnOutAngle = deg2rad(6.0);  // Set the burn out angle (flight-path angle at end of first burn) [rad]
////    const double finalBurnOutMass = 60.7;       // Set the final burn out mass (empty mass + OS mass + excess propellant mass) [kg]
//    const double initialLongitudeDeg = 0.0;     // Set the launch latitude in [deg] (tau)
//    const double initialLatitudeDeg = 0.0;        // Starting latitude [deg] initial condition in (delta)
//    const double HeadingAngle = deg2rad(90.0);  // Set the launch azimuth [rad] (psi)
////    Mars.setUpperAltitudeBound(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the temperature
////    MAV.setUpdatedFinalAltitude(desiredOrbitalAltitude); // Set the upper bound of the altitude [km] for the thrust angles
//    const double initialAltitude = -0.6;          // Starting altitude [km MOLA] initial condition is -0.6 km MOLA
//    const double initialGroundVelocity = 0.0;          // Starting velocity in km/s (is suppose to be 0.0...) 0.00001 default at initial step-size of 0.01 sec



    //////////////////////////////////////////////////////// Test Cases ////////////////////////////////////////////////////////

    /// Comparison?
    const bool comparison = true;

    /// Set initial flight path angle and heading angle
    const double FlightPathAngle = deg2rad(89.0);     // Set flight-path angle in rad --> Default = 90.0 deg
//    const double HeadingAngle = deg2rad(90.0);           // Set heading angle in rad --> Default = 0.0 deg


  /// Initial conditions /// a.k.a. control centre

//    const double initialBurnTime = 68.63; // Burn time from launch till coast
    const double setEndTime = 2000.0;  // Integration end time  // 77 sec for a remainder mass of about 100 kg  // 200 sec for free fall // 2000 for test cases
    const double RKFinitiaterTime = 1.0;    // Time that the RKF integrator is used for TSI initially
    const double EndAltitude = desiredOrbitalAltitude; // Integration end altitude
    const double coastStartTime = initialBurnTime; // Integration coast start time [sec] // test 68.63

//std::cout<<"pi = "<<(4*atan(1))<<std::endl;

    /// TSI settings ///
    const int maxOrder = 20; // Eventually want order 20 (testing is 8)
    /// TSI settings ///

    /// Integration settings ///
    const double chosenLocalErrorTolerance = 1e-8;      // The chosen local error tolerance for TSI
    const double chosenStepSize = 0.01; // The chosen initial step-size for TSI

    std::cout<<"The chosen local error tolerance = "<<chosenLocalErrorTolerance<<std::endl;
    std::cout<<"The chosen initial step-size = "<<chosenStepSize<<std::endl;
    std::cout<<"The chosen end time = "<<setEndTime<<std::endl;
    std::cout<<"The initial Flight-path angle = "<<rad2deg(FlightPathAngle)<<" deg"<<std::endl;
    std::cout<<"The initial Heading angle = "<<rad2deg(HeadingAngle)<<" deg"<<std::endl;

//    const std::string currentIntegrator = "TSI";

    /// Integration settings ///


    // Launch site characteristics


//    const double initialAltitude = -0.6;                 // Starting altitude [km MOLA] initial condition is -0.6 km MOLA
    std::cout<<"The initial altitude = "<<initialAltitude<<std::endl;
//    const double initialLatitudeDeg = 21.0;               // Starting latitude [deg] initial condition is 21 deg (delta)
//    const double initialLongitudeDeg = 74.5;            // Starting longitude [deg] initial condition is 74.5 deg (tau)
//    const double initialGroundVelocity = 0.00001;          // Starting velocity in km/s (is suppose to be 0.0...) 0.00001 default at initial step-size of 0.01 sec
    std::cout<<"The initial ground velocity = "<<initialGroundVelocity<<" km/s"<<std::endl;
    std::cout<<"The initial latitude = "<<initialLatitudeDeg<<" deg"<<std::endl;
    std::cout<<"The initial longitude = "<<initialLongitudeDeg<<" deg"<<std::endl;
    std::cout<<"The order of TSI = "<<maxOrder<<std::endl;



    const double initialLatitude = deg2rad(initialLatitudeDeg);       // Starting latitude [rad]
    const double initialLongitude = deg2rad(initialLongitudeDeg);     // Starting longitude [rad]

    const double initialRadius = bodyReferenceRadius+initialAltitude;               // Starting radius in km

        // Converting the initial spherical position to cartesian position using the standard convertSphericalToCartesian function of Tudat
        // Please note that this function requires the zenith angle as input which is pi/2-latitude!

    Eigen::Vector3d initialCartesianPositionRotationalFrame = Eigen::Vector3d::Zero(3);

      initialCartesianPositionRotationalFrame(0) = initialRadius*cos(initialLatitude)*cos(initialLongitude); // x_R
      initialCartesianPositionRotationalFrame(1) = initialRadius*cos(initialLatitude)*sin(initialLongitude); // y_R
      initialCartesianPositionRotationalFrame(2) = initialRadius*sin(initialLatitude); // z_R

    const Eigen::Vector3d initialCartesianPositionInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocityMars*inertialFrameTime-primeMeridianAngle)*initialCartesianPositionRotationalFrame;


    // Compute initial velocity in y-direction as seen from the launch site in the inertial frame

//    const Eigen::Vector3d initialVelocityLaunchSite = Eigen::Vector3d(0,(rotationalVelocityMars*initialRadius*cos(initialLatitude)),0);

    Eigen::Vector3d initialVelocityLaunchSiteVertical;
    initialVelocityLaunchSiteVertical(0) = initialGroundVelocity*cos(FlightPathAngle)*cos(HeadingAngle);
    initialVelocityLaunchSiteVertical(1) = initialGroundVelocity*cos(FlightPathAngle)*sin(HeadingAngle);
    initialVelocityLaunchSiteVertical(2) = -initialGroundVelocity*sin(FlightPathAngle);

    const Eigen::Vector3d initialVelocityLaunchSite =  tudat::reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(initialLongitude,initialLatitude)*initialVelocityLaunchSiteVertical;

//    std::cout<<"initialVelocityLaunchSite = "<<initialVelocityLaunchSite<<std::endl;

    Eigen::Vector3d initialVelocityDueToMars;

    initialVelocityDueToMars(0) = -rotationalVelocityMars*initialCartesianPositionInertialFrame(1);
    initialVelocityDueToMars(1) = rotationalVelocityMars*initialCartesianPositionInertialFrame(0);
    initialVelocityDueToMars(2) = 0.0;


    /// Debug ///
//    std::cout<<"initialVelocityLaunchSiteVertical = "<<initialVelocityLaunchSiteVertical<<std::endl;
//    std::cout<<"initialVelocityLaunchSite = "<<initialVelocityLaunchSite<<std::endl;
//    std::cout<<"initialVelocityDueToMars = "<<initialVelocityDueToMars<<std::endl;


    /// Debug ///

    const Eigen::Vector3d initialVelocityInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocityMars*inertialFrameTime-primeMeridianAngle)*initialVelocityLaunchSite+initialVelocityDueToMars;

//    const Eigen::Vector3d initialVelocityInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocityMars*inertialFrameTime-primeMeridianAngle+initialLongitude)*initialVelocityLaunchSite;

    /// Setting StateAndTime class using modified vector (Cartesian) ///

    tudat::basic_mathematics::Vector7d aState;

    aState(0) = initialCartesianPositionInertialFrame(0);
    aState(1) = initialCartesianPositionInertialFrame(1);
    aState(2) = initialCartesianPositionInertialFrame(2);
    aState(3) = initialVelocityInertialFrame(0);
    aState(4) = initialVelocityInertialFrame(1);
    aState(5) = initialVelocityInertialFrame(2);
    aState(6) = MAV.MAVmass();  // Mass [kg] from literature study



//    /// Debug ///

//        aState(0) = 0;
//        aState(1) = 0;
//        aState(2) = 1224.18491564675;
//        aState(3) = 0;
//        aState(4) = 0;
//        aState(5) = 0;
//        aState(6) = 227;  // Mass [kg] from literature study

//    /// Debug //




    /// Cartesian state ( in inertial frame) ///

    StateAndTime currentStateAndTime(aState);        // Creating the current state class using the namespace and class directly (Cartesian)

    std::cout<<"aState = "<<aState<<std::endl;

    const tudat::basic_mathematics::Vector7d initialState = aState; // Used for the RKF integrator


    /// Spherical state (in rotating frame) ///

    tudat::basic_mathematics::Vector7d sphericalState;

    sphericalState(0) = initialRadius;   // Radius (r)  [km]
    sphericalState(1) = initialLatitude;   // Latitude (delta)  [rad]
    sphericalState(2) = initialLongitude;   // Longitude (tau)  [rad]
    sphericalState(3) = initialGroundVelocity;   // Ground velocity (V_G)
    sphericalState(4) = FlightPathAngle;   // Flight-path angle (gamma_G)  [rad]
    sphericalState(5) = HeadingAngle;   // Azimuth angle (chi_G)    [rad]
    sphericalState(6) = MAV.MAVmass();    // Mass [kg] from literature study

    StateAndTime currentSphericalStateAndTime(sphericalState);  // Creating the current state class using the namespace and class directly (Spherical)

    std::cout<<"sphericalState = "<<sphericalState<<std::endl;

/////////////////////////////////////////////////////////////////////
////////////////////// Testing the integrators //////////////////////
/////////////////////////////////////////////////////////////////////

/*    /// Debug ///

    const double rotationalVelocityMars = Mars.rotationalVelocity();
    const double primeMeridianAngleMars = Mars.primeMeridianAngle();
    const double inertialFrameTimeMars = Mars.inertialFrameTime();

    const tudat::basic_mathematics::Vector7d currentState = currentStateAndTime.getCurrentState();

    const double xPosition = currentState(0);            // x position coordinate definition
    const double yPosition = currentState(1);            // y position coordinate definition
    const double zPosition = currentState(2);            // z position coordinate definition
    const double xVelocity = currentState(3);            // x velocity coordinate definition
    const double yVelocity = currentState(4);            // y velocity coordinate definition
    const double zVelocity = currentState(5);            // z velocity coordinate definition
    const double massMAV = currentState(6);              // MAV mass definition
    const double currentTime = currentStateAndTime.getCurrentTime();   // current time definition

    const double Radius = sqrt(xPosition*xPosition+yPosition*yPosition+zPosition*zPosition);         // r [km]

    const double inertialVelocity = sqrt(xVelocity*xVelocity+yVelocity*yVelocity+zVelocity*zVelocity);       // V_I [km/s]

    const double inertialLongitude = atan2(yPosition,xPosition);         // lambda [rad]

    const double Latitude = asin(zPosition/Radius);              // delta [rad]

    const double rotationalLongitude = inertialLongitude-rotationalVelocityMars*(inertialFrameTimeMars+currentTime)+primeMeridianAngleMars;  // tau [rad]

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


    const double RadiusChange = (xPosition*xVelocity+yPosition*yVelocity+zPosition*zVelocity)/(Radius);      // radial velocity [km/s]

    double LatitudeChange_; // delta_dot [rad/s] (placeholder)

    if ((Radius*Radius*sqrt(1-(zPosition/Radius)*(zPosition/Radius))) == 0){
        LatitudeChange_ = 0;
    }
    else{
        LatitudeChange_ = (Radius*zVelocity-zPosition*RadiusChange)/(Radius*Radius*sqrt(1-(zPosition/Radius)*(zPosition/Radius)));
    };

    const double LatitudeChange = LatitudeChange_;  // delta_dot [rad/s] (actual parameter)

    const double localMarsRotationalVelocity = rotationalVelocityMars*Radius*cos(Latitude);  // V_M [km/s]

    const double inertialFlightPathAngle = asin(RadiusChange/inertialVelocity);          // gamma_I [rad]

    const double inertialAzimuth = atan2((inertialLongitudeChange*cos(Latitude)),LatitudeChange);    // chi_I [rad]

    const double rotationalVelocityMAV = sqrt(localMarsRotationalVelocity*localMarsRotationalVelocity+inertialVelocity*inertialVelocity-2*localMarsRotationalVelocity*inertialVelocity*cos(inertialFlightPathAngle)*sin(inertialAzimuth));  // V_R [km/s]

    double rotationalFlightPathAngle_; // gamma_R [rad]  (placeholder)

    if (rotationalVelocityMAV == 0){       // Setting the initial flight path angle in the rotational frame to 90 deg (or pi/s)

        rotationalFlightPathAngle_ = tudat::mathematical_constants::LONG_PI/2;
    }
    else {
        rotationalFlightPathAngle_ = asin(RadiusChange/rotationalVelocityMAV);
    };

    const double rotationalFlightPathAngle = rotationalFlightPathAngle_;    // gamma_R [rad] (actual parameter)


    const double rotationalAzimuth = atan2((rotationalLongitudeChange*cos(Latitude)),LatitudeChange);    // chi_R [rad]

    // Check output
        std::cout<<"Radius = "<<Radius<<std::endl;
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
        std::cout<<"rotationalVelocityMAV = "<<rotationalVelocityMAV<<std::endl;
        std::cout<<"rotationalFlightPathAngle = "<<rotationalFlightPathAngle<<std::endl;
        std::cout<<"rotationalAzimuth = "<<rotationalAzimuth<<std::endl;



        /// Testing the local air temperature function ///

            const double currentAltitude = Radius-Mars.bodyReferenceRadius();

            const double currentTemperature = air_temperature::airTemperature(Mars.temperaturePolyCoefficients(), Mars.temperatureAltitudeRanges(),currentAltitude);

           // Check output
            std::cout<<"currentAltitude = "<<currentAltitude<<std::endl;
            std::cout<<"currentTemperature = "<<currentTemperature<<std::endl;
            std::cout<<"Radius = "<<Radius<<std::endl;
//            std::cout<<"R_MOLA = "<<Mars.bodyReferenceRadius()<<std::endl;
//            std::cout<<"Radius-3395.4 = "<<Radius-3395.4<<std::endl;
//            std::cout<<"R_MOLA-3396 = "<<Mars.bodyReferenceRadius()-3396<<std::endl;


        /// Testing the local air density function ///

            const double currentDensity= air_density::airDensity(Mars.densityPolyCoefficients(),  currentAltitude);

            std::cout<<"The current air density = "<<currentDensity<<std::endl;

        /// Testing the ascentDragForce function ///



            const double currentDrag = Drag::ascentDragForce(rotationalVelocityMAV,currentTemperature,Mars.adiabeticIndex(),Mars.specificGasConstant(),
                                                             MAV.dragCoefficientPolyCoefficients(),MAV.dragCoefficientMachRanges(),
                                                             MAV.referenceArea(), currentDensity);

            std::cout<<"currentDrag = "<<currentDrag<<std::endl;

    /// Debug ///
//*/


//////////////////////////////////////////////////////////////////////////////////
////////////////////// Testing the Taylor series integrator //////////////////////
//////////////////////////////////////////////////////////////////////////////////

/// Setting the data collection file for TSI and inserting the first values ///

    // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
    const std::string outputDirectoryTSI = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/01.integrationResults/TSI/";


    // Set output format for matrix output.
    Eigen::IOFormat csvFormatTSI( 15, 0, ", ", "\n" );

    // Set absolute path to file containing the Taylor Series Coefficients.
    std::string dataAbsolutePathTSI = outputDirectoryTSI + "test6FullIntegrationTSIstateAndTime.csv";
    std::string dataAbsolutePathTSICart = outputDirectoryTSI + "test6FullIntegrationTSIstateAndTime.csv";

    // Create a row vector for the storing of the data
    Eigen::MatrixXd outputVectorTSI = Eigen::MatrixXd::Zero(1,8); // Create a row vector for the storing of the data
    Eigen::MatrixXd outputVectorTSICartesian = Eigen::MatrixXd::Zero(1,8); // Create a row vector for the storing of the cartesian data

    // Getting the initial conditions for storage
    const tudat::basic_mathematics::Vector7d initialStateTSI = currentSphericalStateAndTime.getCurrentState();
    const tudat::basic_mathematics::Vector7d initialStateTSICartesian = currentStateAndTime.getCurrentState();

    // Filling the output vector
    outputVectorTSI(0,0) = currentSphericalStateAndTime.getCurrentTime();   // Storing the initial time
    outputVectorTSI(0,1) = initialStateTSI(0);   // Storing the initial radius
    outputVectorTSI(0,2) = initialStateTSI(1);   // Storing the initial latitude
    outputVectorTSI(0,3) = initialStateTSI(2);   // Storing the initial longitude
    outputVectorTSI(0,4) = initialStateTSI(3);   // Storing the initial ground velocity
    outputVectorTSI(0,5) = initialStateTSI(4);   // Storing the initial flight-path angle
    outputVectorTSI(0,6) = initialStateTSI(5);   // Storing the initial azimuth angle
    outputVectorTSI(0,7) = initialStateTSI(6);   // Storing the initial MAV mass

    outputVectorTSICartesian(0,0) = currentStateAndTime.getCurrentTime();   // Storing the initial time
    outputVectorTSICartesian(0,1) = initialStateTSICartesian(0);   // Storing the initial x position
    outputVectorTSICartesian(0,2) = initialStateTSICartesian(1);   // Storing the initial y position
    outputVectorTSICartesian(0,3) = initialStateTSICartesian(2);   // Storing the initial z position
    outputVectorTSICartesian(0,4) = initialStateTSICartesian(3);   // Storing the initial x velocity
    outputVectorTSICartesian(0,5) = initialStateTSICartesian(4);   // Storing the initial y velocity
    outputVectorTSICartesian(0,6) = initialStateTSICartesian(5);   // Storing the initial z velocity
    outputVectorTSICartesian(0,7) = initialStateTSICartesian(6);   // Storing the initial MAV mass

    // Storing the data

    std::ifstream ifileTSI(dataAbsolutePathTSI.c_str()); // Check it as an input file

    bool fexistsTSI = false;   // Set the default to "It does not exist"

    if (ifileTSI){         // Attempt to open the file


       fexistsTSI = true;      // If the file can be opened it must exist

       ifileTSI.close();   // Close the file

    }


    // If so: error and create temporary file, if not: create new file and put data in

    if (fexistsTSI == true){


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

        std::string ComputerTimeString = currentYear + "-" + currentMonth + "-" + currentDay + "_" + currentHour + ":" + currentMin + ":" + currentSec;  // Convert to string and store

        std::string newFileName = "backupSphericalTSIFileAtDateAndTime_" + ComputerTimeString + ".csv";
        std::string newFileNameCart = "backupSpherical(Cart)TSIFileAtDateAndTime_" + ComputerTimeString + ".csv";

    std::cerr<<"The file name that you have chosen already exists, a new file with name "<<newFileName<<" will be created to store the data for now"<<std::endl;
    std::cerr<<"The file name that you have chosen already exists, a new file with name "<<newFileNameCart<<" will be created to store the data for now"<<std::endl;

    // Set new absolute path to file containing the data.
    dataAbsolutePathTSI = outputDirectoryTSI + newFileName;
    dataAbsolutePathTSICart = outputDirectoryTSI + newFileNameCart;

    // Export the data.
    std::ofstream exportFile1( dataAbsolutePathTSI.c_str( ) ); // Make the new file
    std::cout<<"New file called "<<dataAbsolutePathTSI<<" has been created"<<std::endl;
    exportFile1 << outputVectorTSI.format( csvFormatTSI );          // Store the new values
    exportFile1.close( );   // Close the file

    std::ofstream exportFile1Cart( dataAbsolutePathTSICart.c_str( ) ); // Make the new file
    std::cout<<"New file called "<<dataAbsolutePathTSICart<<" has been created"<<std::endl;
    exportFile1Cart << outputVectorTSICartesian.format( csvFormatTSI );          // Store the new values
    exportFile1Cart.close( );   // Close the file


}
        else{

        // Export the data.
        std::ofstream exportFile1( dataAbsolutePathTSI.c_str( ) ); // Make the new file
        std::cout<<"New file called "<<dataAbsolutePathTSI<<" has been created"<<std::endl;
        exportFile1 << outputVectorTSI.format( csvFormatTSI );          // Store the new values
        exportFile1.close( );   // Close the file

        std::ofstream exportFile1Cart( dataAbsolutePathTSICart.c_str( ) ); // Make the new file
        std::cout<<"New file called "<<dataAbsolutePathTSICart<<" has been created"<<std::endl;
        exportFile1Cart << outputVectorTSICartesian.format( csvFormatTSI );          // Store the new values
        exportFile1Cart.close( );   // Close the file
    };


    // Define storing matrix for the intermediate values
    Eigen::MatrixXd dataStoringMatrixTSI(1,8); // The size of this matrix will change in the do-loop
    Eigen::MatrixXd dataStoringMatrixTSICartesian(1,8); // The size of this matrix will change in the do-loop
    tudat::basic_mathematics::Vector7d stateAtPoint2SecTSI; // Storing the 0.2 seconds value specifically for comparison


    ///// First steps by RKF integrator /////


    /// Testing with the state derivative function class


    // Full complete test

    ascentStateDerivativeFunctionClass stateDerivativeFunctionClassTSI(Mars,MAV);     // Initialize the class

    // Set the initial values for the flight-path angle and heading angle
    stateDerivativeFunctionClassTSI.setFlightPathAngleAndHeadingAngle(FlightPathAngle,HeadingAngle);

    // Just using boost::bind

    boost::function< tudat::basic_mathematics::Vector7d( const double, const tudat::basic_mathematics::Vector7d ) > stateDerivativeFunctionTSI // Using boost::function to create the passing function
            =boost::bind( &ascentStateDerivativeFunctionClass::ascentStateDerivativeFunction, &stateDerivativeFunctionClassTSI, _1, _2);       // Then using boost::bind to bind the function as per class stateDerivativeFunctionClass which requires two inputs:
                                                                                                                                            // _1 which is the current time and _2 which is the current state


 ///// Testing the implementation in the integrator ///

    // Initial conditions.
//                    const double initialTime = currentStateAndTime.getCurrentTime();                            // Time.
    const double initialTimeTSI = 0;                            // Time. set for verification
//    Eigen::VectorXd initialState = currentStateAndTime.getCurrentState(); // State: start with zero velocity at the origin.

    const double endTimeTSIRKF = setEndTime;     // Using the same initial step-size as defined for TSI RKF

    // Step-size settings.
    // The minimum and maximum step-size are set such that the input data is fully accepted by the
    // integrator, to determine the steps to be taken.
    const double zeroMinimumStepSizeTSI = std::numeric_limits< double >::epsilon( );
    const double infiniteMaximumStepSizeTSI = std::numeric_limits< double >::infinity( );
    double stepSizeRKFTSI = chosenStepSize;          // Using the same initial step-size as defined for TSI

    // Tolerances.
    const double relativeToleranceTSI = 1e-15; /*chosenLocalErrorTolerance;     //*/
    const double absoluteToleranceTSI = 1e-15; /*chosenLocalErrorTolerance;     //*/





////////////////////////// RKF 7(8) integrator is used in this case///////////////////////////////////////////////////////////////////////////////////////////////////////



    /// Runge-Kutta-Fehlberg 7(8) integrator.
       tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integratorTSI(
                   tudat::numerical_integrators::RungeKuttaCoefficients::get(
                       tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78),
                   stateDerivativeFunctionTSI, initialTimeTSI, initialState, zeroMinimumStepSizeTSI,
                   infiniteMaximumStepSizeTSI, relativeToleranceTSI, absoluteToleranceTSI );



//       /// Runge-Kutta-Fehlberg 4(5) integrator.
//           tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integratorTSI(
//                       tudat::numerical_integrators::RungeKuttaCoefficients::get(
//                           tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg45),
//                       stateDerivativeFunctionTSI, initialTimeTSI, initialState, zeroMinimumStepSizeTSI,
//                       infiniteMaximumStepSizeTSI, relativeToleranceTSI, absoluteToleranceTSI );

//                          /// Dormand-Prince 8(7) integrator.
//           tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integratorTSI(
//                       tudat::numerical_integrators::RungeKuttaCoefficients::get(
//                           tudat::numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince),
//                       stateDerivativeFunctionTSI, initialTimeTSI, initialState, zeroMinimumStepSizeTSI,
//                       infiniteMaximumStepSizeTSI, relativeToleranceTSI, absoluteToleranceTSI );




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Set initial running time. This is updated after each step that the numerical integrator takes.
    double runningTimeTSI = 0.0;
    int countRKFTSI = 0;


    // Define storing matrix for the intermediate values
//    Eigen::MatrixXd dataStoringMatrix(1,8); // The size of this matrix will change in the do-loop


    do
    {
        // Make sure the integrator does not integrate beyond the end time.



        if ( std::fabs( endTimeTSIRKF - runningTimeTSI )
             <= std::fabs( stepSizeRKFTSI ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) )
        {
            stepSizeRKFTSI = endTimeTSIRKF - integratorTSI.getCurrentIndependentVariable( );
        }

//        double prevStepSizeTSI = stepSizeRKFTSI;

//                         std::cout<<"The current stepSize is "<<prevStepSizeTSI<<" s"<<std::endl;

        // Perform a single integration step. Then update the step-size and running time.
        integratorTSI.performIntegrationStep( stepSizeRKFTSI );
        stepSizeRKFTSI = integratorTSI.getNextStepSize( );
        runningTimeTSI = integratorTSI.getCurrentIndependentVariable( );

        Eigen::VectorXd currentStateRKFTSI = integratorTSI.getCurrentState();

//        std::cout<<"currentStateRKFTSI = "<<currentStateRKFTSI<<std::endl;
//        std::cout<<"r = "<<sqrt(currentStateRKFTSI(0)*currentStateRKFTSI(0)+currentStateRKFTSI(1)*currentStateRKFTSI(1)+currentStateRKFTSI(2)*currentStateRKFTSI(2))<<std::endl;


        /// Get the spherical state from the initial few seconds of RKF

        const double angleItoR_1 = rotationalVelocityMars*(inertialFrameTime+runningTimeTSI)-primeMeridianAngle;





        // Position

        Eigen::Vector3d CartesianPositionInertialFrameRKF = Eigen::Vector3d::Zero(1,3); // Define the cartesian position inertial frame
        CartesianPositionInertialFrameRKF(0) = currentStateRKFTSI(0); // x
        CartesianPositionInertialFrameRKF(1) = currentStateRKFTSI(1); // y
        CartesianPositionInertialFrameRKF(2) = currentStateRKFTSI(2); // z


        const Eigen::Vector3d CartesianPositionRotationalFrameRKF = tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix(angleItoR_1)*CartesianPositionInertialFrameRKF;



        Eigen::MatrixXd updatedStateAndTimeVectorRKFTSI = Eigen::MatrixXd::Zero(1,8); // Creating a temporary vector to hold the converted state

        updatedStateAndTimeVectorRKFTSI(0) = sqrt(currentStateRKFTSI(0)*currentStateRKFTSI(0)+currentStateRKFTSI(1)*currentStateRKFTSI(1)+currentStateRKFTSI(2)*currentStateRKFTSI(2)); // Radius

        updatedStateAndTimeVectorRKFTSI(1) = asin(currentStateRKFTSI(2)/updatedStateAndTimeVectorRKFTSI(0)) ; // latitude (delta)

        updatedStateAndTimeVectorRKFTSI(2) = atan2(CartesianPositionRotationalFrameRKF(1),CartesianPositionRotationalFrameRKF(0)); // longitude (tau)




        /// Debug ///

//        std::cout<<"/// Debug ///"<<std::endl;
//        std::cout<<"CartesianPositionInertialFrameRKF = "<<CartesianPositionInertialFrameRKF<<std::endl;
//        std::cout<<"CartesianPositionRotationalFrameRKF = "<<CartesianPositionRotationalFrameRKF<<std::endl;
//        std::cout<<"angleItoR_1 = "<<angleItoR_1<<std::endl;
//        std::cout<<"atan2(currentStateRKFTSI(1),currentStateRKFTSI(0))-updatedStateAndTimeVectorRKFTSI(2) = "<<atan2(currentStateRKFTSI(1),currentStateRKFTSI(0))-updatedStateAndTimeVectorRKFTSI(2)<<std::endl;
//        std::cout<<"atan2(currentStateRKFTSI(1),currentStateRKFTSI(0)) = "<<atan2(currentStateRKFTSI(1),currentStateRKFTSI(0))<<std::endl;
////        std::cout<<"-rotationalVelocityMars*(runningTimeTSI+prevStepSizeTSI) = "<<-rotationalVelocityMars*(runningTimeTSI+prevStepSizeTSI)<<std::endl;
////        std::cout<<"primeMeridianAngle = "<<primeMeridianAngle<<std::endl;



//        std::cout<<"/// Debug ///"<<std::endl;
//        std::cout<<" "<<std::endl;

        /// Debug ///



        // Velocity

        Eigen::Vector3d CartesianVelocityInertialFrameRKF = Eigen::Vector3d::Zero(1,3); // Define the cartesian velocity inertial frame
        CartesianVelocityInertialFrameRKF(0) = currentStateRKFTSI(3); // Vx
        CartesianVelocityInertialFrameRKF(1) = currentStateRKFTSI(4); // Vy
        CartesianVelocityInertialFrameRKF(2) = currentStateRKFTSI(5); // Vz

        Eigen::Vector3d rotatingPlanetCorrection = Eigen::Vector3d::Zero(3);

        rotatingPlanetCorrection(0) = -rotationalVelocityMars*currentStateRKFTSI(1);
        rotatingPlanetCorrection(1) = rotationalVelocityMars*currentStateRKFTSI(0);
        rotatingPlanetCorrection(2) = 0.0;


        const Eigen::Vector3d CartesianVelocityRotationalFrameRKF = tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix(angleItoR_1)*(CartesianVelocityInertialFrameRKF-rotatingPlanetCorrection); // Cartesian velocity in the rotational frame

        const Eigen::Vector3d CartesianVelocityVerticalFrameRKF = tudat::reference_frames::getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(updatedStateAndTimeVectorRKFTSI(2),updatedStateAndTimeVectorRKFTSI(1))*CartesianVelocityRotationalFrameRKF; // Cartesian velocity in the local vertical frame


        updatedStateAndTimeVectorRKFTSI(3) = sqrt((currentStateRKFTSI(3)+rotationalVelocityMars*currentStateRKFTSI(1))*(currentStateRKFTSI(3)+rotationalVelocityMars*currentStateRKFTSI(1))+(currentStateRKFTSI(4)-rotationalVelocityMars*currentStateRKFTSI(0))*(currentStateRKFTSI(4)-rotationalVelocityMars*currentStateRKFTSI(0))+currentStateRKFTSI(5)*currentStateRKFTSI(5)); // Ground velocity scalar


//        const double V_x_V_RKFTSI = -(currentStateRKFTSI(3)+rotationalVelocityMars*currentStateRKFTSI(1))*sin(updatedStateAndTimeVectorRKFTSI(1))*cos(atan2(currentStateRKFTSI(1),currentStateRKFTSI(0)))-
//                (currentStateRKFTSI(4)-rotationalVelocityMars*currentStateRKFTSI(0))*sin(updatedStateAndTimeVectorRKFTSI(1))*sin(atan2(currentStateRKFTSI(1),currentStateRKFTSI(0)))+
//                currentStateRKFTSI(5)*cos(updatedStateAndTimeVectorRKFTSI(1)); // Vxv

//        const double V_y_V_RKFTSI = (currentStateRKFTSI(4)-rotationalVelocityMars*currentStateRKFTSI(0))*cos(atan2(currentStateRKFTSI(1),currentStateRKFTSI(0)))
//                -(currentStateRKFTSI(3)+rotationalVelocityMars*currentStateRKFTSI(1))*sin(atan2(currentStateRKFTSI(1),currentStateRKFTSI(0))); // Vyv

//        const double V_z_V_RKFTSI = -(currentStateRKFTSI(3)+rotationalVelocityMars*currentStateRKFTSI(1))*cos(updatedStateAndTimeVectorRKFTSI(1))*cos(atan2(currentStateRKFTSI(1),currentStateRKFTSI(0)))-
//                (currentStateRKFTSI(4)-rotationalVelocityMars*currentStateRKFTSI(0))*cos(updatedStateAndTimeVectorRKFTSI(1))*sin(atan2(currentStateRKFTSI(1),currentStateRKFTSI(0)))-
//                currentStateRKFTSI(5)*sin(updatedStateAndTimeVectorRKFTSI(1)); // Vzv


        /// Debug ///
//        std::cout<<"/// Debug ///"<<std::endl;
//        std::cout<<"CartesianVelocityVerticalFrameRKF = "<<CartesianVelocityVerticalFrameRKF<<std::endl;
//        std::cout<<"V_x_V_RKFTSI = "<<V_x_V_RKFTSI<<std::endl;
//        std::cout<<"V_y_V_RKFTSI = "<<V_y_V_RKFTSI<<std::endl;
//        std::cout<<"V_z_V_RKFTSI = "<<V_z_V_RKFTSI<<std::endl;
//        std::cout<<"Lambda = "<<updatedStateAndTimeVectorRKFTSI(2)+rotationalVelocityMars*((runningTimeTSI+prevStepSizeTSI))<<std::endl;
//        std::cout<<"Lambda (from equation) = "<<
//        std::cout<<"/// Debug ///"<<std::endl;
//        std::cout<<" "<<std::endl;

        /// Debug ///


//        updatedStateAndTimeVectorRKFTSI(4) = -asin(V_z_V_RKFTSI/updatedStateAndTimeVectorRKFTSI(3)); // flight-path angle

//        std::cout<<"flight-path angle = "<<updatedStateAndTimeVectorRKFTSI(4)<<std::endl;


//        updatedStateAndTimeVectorRKFTSI(5) = atan2(V_y_V_RKFTSI,V_x_V_RKFTSI); // azimuth angle

//        std::cout<<"azimuth angle = "<<updatedStateAndTimeVectorRKFTSI(5)<<std::endl;

        updatedStateAndTimeVectorRKFTSI(4) = -asin(CartesianVelocityVerticalFrameRKF(2)/updatedStateAndTimeVectorRKFTSI(3)); // flight-path angle

//        std::cout<<"flight-path angle = "<<updatedStateAndTimeVectorRKFTSI(4)<<std::endl;


        updatedStateAndTimeVectorRKFTSI(5) = atan2(CartesianVelocityVerticalFrameRKF(1),CartesianVelocityVerticalFrameRKF(0)); // azimuth angle

//        std::cout<<"azimuth angle = "<<updatedStateAndTimeVectorRKFTSI(5)<<std::endl;

//        std::cout<<" "<<std::endl;


        updatedStateAndTimeVectorRKFTSI(6) = currentStateRKFTSI(6); // MAV mass

//        std::cout<<"updatedStateAndTimeVectorRKFTSI = "<<updatedStateAndTimeVectorRKFTSI<<std::endl;

        outputVectorTSI = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.

        // Filling the output vector
        outputVectorTSI(0,0) = runningTimeTSI;   // Storing the updated time
        outputVectorTSI(0,1) = updatedStateAndTimeVectorRKFTSI(0);   // Storing the updated radius
        outputVectorTSI(0,2) = updatedStateAndTimeVectorRKFTSI(1);   // Storing the updated latitude
        outputVectorTSI(0,3) = updatedStateAndTimeVectorRKFTSI(2);   // Storing the updated longitude
        outputVectorTSI(0,4) = updatedStateAndTimeVectorRKFTSI(3);   // Storing the updated ground velocity
        outputVectorTSI(0,5) = updatedStateAndTimeVectorRKFTSI(4);   // Storing the updated flight-path angle
        outputVectorTSI(0,6) = updatedStateAndTimeVectorRKFTSI(5);   // Storing the updated azimuth angle
        outputVectorTSI(0,7) = updatedStateAndTimeVectorRKFTSI(6);   // Storing the updated MAV mass

        outputVectorTSICartesian = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.

        // Filling the output vector
        outputVectorTSICartesian(0,0) = runningTimeTSI;   // Storing the updated time
        outputVectorTSICartesian(0,1) = currentStateRKFTSI(0);   // Storing the updated x-position
        outputVectorTSICartesian(0,2) = currentStateRKFTSI(1);   // Storing the updated y-position
        outputVectorTSICartesian(0,3) = currentStateRKFTSI(2);   // Storing the updated z-position
        outputVectorTSICartesian(0,4) = currentStateRKFTSI(3);   // Storing the updated x-velocity
        outputVectorTSICartesian(0,5) = currentStateRKFTSI(4);   // Storing the updated y-velocity
        outputVectorTSICartesian(0,6) = currentStateRKFTSI(5);   // Storing the updated z-velocity
        outputVectorTSICartesian(0,7) = currentStateRKFTSI(6);   // Storing the updated MAV mass

        // Store the new values in the data storage matrix

        if (countRKFTSI == 0){

          dataStoringMatrixTSI.row(countRKFTSI) = outputVectorTSI.row(0); // Filling the matrix
          dataStoringMatrixTSICartesian.row(countRKFTSI) = outputVectorTSICartesian.row(0); // Filling the matrix

        }
        else{
            dataStoringMatrixTSI.conservativeResize(countRKFTSI+1,8); // Making the matrix bigger in order to store more values
            dataStoringMatrixTSICartesian.conservativeResize(countRKFTSI+1,8); // Making the matrix bigger in order to store more values

            dataStoringMatrixTSI.row(countRKFTSI) = outputVectorTSI.row(0); // Filling the matrix
            dataStoringMatrixTSICartesian.row(countRKFTSI) = outputVectorTSICartesian.row(0); // Filling the matrix
        }

        // Updating the current state and time class!!!
        tudat::basic_mathematics::Vector7d currentSphericalStateVector; // Create the current state and time vector

        // Fill the curent state and time vector
        currentSphericalStateVector(0) = updatedStateAndTimeVectorRKFTSI(0);  // Updated radius
        currentSphericalStateVector(1) = updatedStateAndTimeVectorRKFTSI(1);  // Updated latitude
        currentSphericalStateVector(2) = updatedStateAndTimeVectorRKFTSI(2);  // Updated longitude
        currentSphericalStateVector(3) = updatedStateAndTimeVectorRKFTSI(3);  // Updated ground velocity
        currentSphericalStateVector(4) = updatedStateAndTimeVectorRKFTSI(4);  // Updated flight-path angle
        currentSphericalStateVector(5) = updatedStateAndTimeVectorRKFTSI(5);  // Updated azimuth angle
        currentSphericalStateVector(6) = updatedStateAndTimeVectorRKFTSI(6);  // Updated MAV mass

//        runningTimeTSI = updatedStateAndTimeVector(7);             // Updated time

//        std::cout<<"currentSphericalStateVector = "<<currentSphericalStateVector<<std::endl;




        currentSphericalStateAndTime.setCurrentStateAndTime(currentSphericalStateVector,runningTimeTSI); // Update the current state and time class!

        countRKFTSI++;

//        }while( !( endTimeTSIRKF - runningTimeTSI <= std::numeric_limits< double >::epsilon( ) ) );
        }while( ( runningTimeTSI <= RKFinitiaterTime) );

        ///// First steps by RKF integrator /////

    // Input for RKF integrator (after TSI)
  const double finalTimeFirstSecond = runningTimeTSI;
   const tudat::basic_mathematics::Vector7d finalStateFirstSecond = integratorTSI.getCurrentState();
   const double finalStepSizeFirstSecond = stepSizeRKFTSI;



//        std::cout<<"dataStoringMatrixTSI = "<<dataStoringMatrixTSI<<std::endl;


/// Defining the order and initializing the StepSize class ///

std::cout<<"////////////////////////////////////////////////////////////////// Start of TSI //////////////////////////////////////////////////////////////////"<<std::endl;


        StepSize stepSize; // Initializing the stepSize class. THIS SHOULD BE DONE BEFORE THE START OF THE INTEGRATION!!!!!

        stepSize.setLocalErrorTolerance(chosenLocalErrorTolerance);  // Setting the local error tolerance to the lowest possible value in order to compare to RKF7(8) and the others
//        stepSize.setCurrentStepSize(chosenStepSize); // Setting the step-size to the chosen step-size
        stepSize.setCurrentStepSize(stepSizeRKFTSI); // Setting the step-size to the chosen step-size


//        /// Debug ///

//        std::cout<<"StepSize has been initialized"<<std::endl;
//        stepSize.setCurrentStepSize(1);
//        std::cout<<"The current step-size = "<<stepSize.getCurrentStepSize()<<std::endl;
//        stepSize.setCurrentStepSize(5);
//        std::cout<<"The current step-size = "<<stepSize.getCurrentStepSize()<<std::endl;

//        /// Debug ///

//        const double currentStepSize = stepSize.getCurrentStepSize();

//        std::cout<<"The current stepSize = "<<currentStepSize<<std::endl;


/// Performing the actual TSI integration ///

        // Determine the CPU time

        const double initialCPUTime = clock();
        // Determine the CPU time


        // Set the end time
        const double endTimeTSI = setEndTime; // sec

//        std::cout<<"It works till here 1"<<std::endl;


/// The integeration do-loop ///

        // Set initial running time. This is updated after each step that the numerical integrator takes.
//     double runningTimeTSI = 0.0;
     int countTSI = 0;
     bool coast = false;    // to determine if there is going to be a coasting phase

     if (coastStartTime<endTimeTSI){
         coast = true;
     }

//     std::cout<<"coast = "<<coast<<std::endl;

    do
    {

//     std::cout<<"It works till here 2"<<std::endl;

//     for (int i = 0; i<4; i++){
         /// Debug ///
//    std::cout<<"The current step-size = "<<stepSize.getCurrentStepSize()<<std::endl;
//    std::cout<<"The current runningTime = "<<runningTimeTSI<<std::endl;
//    std::cout<<"std::fabs(endTime-runningTime) = "<<std::fabs(endTimeTSI-runningTimeTSI)<<std::endl;
//    std::cout<<"std::fabs( stepSize.getCurrentStepSize() ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) = "<<std::fabs( stepSize.getCurrentStepSize() ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) )<<std::endl;
         /// Debug ///

//    stepSize.setCurrentStepSize(25); // Specifying a constant step-size for verification
//         std::cout<<"runningTimeTSI = "<<runningTimeTSI<<std::endl;

         if (coast == true && runningTimeTSI < coastStartTime){
             if ( std::fabs( coastStartTime - runningTimeTSI )
                                      <= std::fabs( stepSize.getCurrentStepSize() ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) )
                                 {
//                                        std::cout<<"Current stepSize = "<<stepSize.getCurrentStepSize()<<std::endl;
                                     stepSize.setCurrentStepSize(coastStartTime - runningTimeTSI);
//                                     coast = false;
//                                     std::cout<<"This should only happen once! And the stepSize = "<<stepSize.getCurrentStepSize()<<std::endl;

                                 }
//             std::cout<<"This should be happening all the time!"<<std::endl;

         } // if coast
         else {

//              std::cout<<"Current stepSize = "<<stepSize.getCurrentStepSize()<<std::endl;
//              std::cout<<"endTimeTSI-runningTimeTSI = "<<fabs(endTimeTSI-runningTimeTSI)<<std::endl;

    if ( std::fabs( endTimeTSI - runningTimeTSI )
                             <= std::fabs( stepSize.getCurrentStepSize() ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) )
                        {
//                            std::cout<<"Why does it not go here?"<<std::endl;
                            stepSize.setCurrentStepSize(endTimeTSI - runningTimeTSI);
                        }

} // if not coast

//          std::cout<<"Current stepSize = "<<stepSize.getCurrentStepSize()<<std::endl;
//          std::cout<<"Updated time = "<<runningTimeTSI+stepSize.getCurrentStepSize()<<std::endl;

        Eigen::VectorXd updatedStateAndTimeVector = performTaylorSeriesIntegrationStep(Mars, MAV, currentSphericalStateAndTime, stepSize, maxOrder, FlightPathAngle, HeadingAngle); /// The actual integration step
        // This function has the output: updated position, updated velocity, updated mass and updated time

        if (updatedStateAndTimeVector(0)-bodyReferenceRadius>=320.0){ // Accounting for the case where the altitude goes above 320.0 km MOLA, then assume drag == 0.0 N (density curve was fitted till 320.0 km MOLA)
            MAV.setReferenceArea(0.0);
        }

//        std::cout<<"updatedSphericalStateAndTimeVector = "<<updatedStateAndTimeVector<<std::endl;


        // Check to see if the class has been updated from within the TSI function
//        const double nextStepSize = stepSize.getCurrentStepSize();

//        std::cout<<"The next stepSize = "<<nextStepSize<<std::endl;

/// Storing the values ///

        outputVectorTSI = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.

        // Filling the output vector
        outputVectorTSI(0,0) = updatedStateAndTimeVector(7);   // Storing the updated time
        outputVectorTSI(0,1) = updatedStateAndTimeVector(0);   // Storing the updated radius
        outputVectorTSI(0,2) = updatedStateAndTimeVector(1);   // Storing the updated latitude
        outputVectorTSI(0,3) = updatedStateAndTimeVector(2);   // Storing the updated longitude
        outputVectorTSI(0,4) = updatedStateAndTimeVector(3);   // Storing the updated ground velocity
        outputVectorTSI(0,5) = updatedStateAndTimeVector(4);   // Storing the updated flight-path angle
        outputVectorTSI(0,6) = updatedStateAndTimeVector(5);   // Storing the updated azimuth angle
        outputVectorTSI(0,7) = updatedStateAndTimeVector(6);   // Storing the updated MAV mass




        // Store the new values in the data storage matrix

        if (countTSI+countRKFTSI == 0){

          dataStoringMatrixTSI.row(countTSI+countRKFTSI) = outputVectorTSI.row(0); // Filling the matrix
        }
        else{
            dataStoringMatrixTSI.conservativeResize(countTSI+countRKFTSI+1,8); // Making the matrix bigger in order to store more values

            dataStoringMatrixTSI.row(countTSI+countRKFTSI) = outputVectorTSI.row(0); // Filling the matrix
        }

        // Updating the current state and time class!!!
        tudat::basic_mathematics::Vector7d currentSphericalStateVector; // Create the current state and time vector

        // Fill the curent state and time vector
        currentSphericalStateVector(0) = updatedStateAndTimeVector(0);  // Updated radius
        currentSphericalStateVector(1) = updatedStateAndTimeVector(1);  // Updated latitude
        currentSphericalStateVector(2) = updatedStateAndTimeVector(2);  // Updated longitude
        currentSphericalStateVector(3) = updatedStateAndTimeVector(3);  // Updated ground velocity
        currentSphericalStateVector(4) = updatedStateAndTimeVector(4);  // Updated flight-path angle
        currentSphericalStateVector(5) = updatedStateAndTimeVector(5);  // Updated azimuth angle
        currentSphericalStateVector(6) = updatedStateAndTimeVector(6);  // Updated MAV mass

        runningTimeTSI = updatedStateAndTimeVector(7);             // Updated time




        currentSphericalStateAndTime.setCurrentStateAndTime(currentSphericalStateVector,runningTimeTSI); // Update the current state and time class!

        /// Compute the inertial cartesian state ///

        Eigen::Vector3d positionRotatingFrame;

        positionRotatingFrame(0) = currentSphericalStateVector(0)*cos(currentSphericalStateVector(1))*cos(currentSphericalStateVector(2));    // x-position
        positionRotatingFrame(1) = currentSphericalStateVector(0)*cos(currentSphericalStateVector(1))*sin(currentSphericalStateVector(2));    // y-position
        positionRotatingFrame(2) = currentSphericalStateVector(0)*sin(currentSphericalStateVector(1));    // z-position

        Eigen::Vector3d velocityVerticalFrame;

        velocityVerticalFrame(0) = currentSphericalStateVector(3)*cos(currentSphericalStateVector(4))*cos(currentSphericalStateVector(5));    // x-velocity
        velocityVerticalFrame(1) = currentSphericalStateVector(3)*cos(currentSphericalStateVector(4))*sin(currentSphericalStateVector(5));    // y-velocity
        velocityVerticalFrame(2) = -currentSphericalStateVector(3)*sin(currentSphericalStateVector(4));    // z-velocity

        // Perform the reference frame transformations

        const double angleItoR = rotationalVelocityMars*(inertialFrameTime+runningTimeTSI)-primeMeridianAngle;

        // Position from R-frame to I-frame
        const Eigen::Vector3d positionInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(angleItoR)*positionRotatingFrame;


        // Velocity from V-frame to R-frame
        const Eigen::Vector3d velocityRotationalFrame = tudat::reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(currentSphericalStateVector(2),currentSphericalStateVector(1))*velocityVerticalFrame;

        // Vector to account for the rotation of Mars
        Eigen::Vector3d velocityDueToMars;

        velocityDueToMars(0) = -rotationalVelocityMars*positionInertialFrame(1);
        velocityDueToMars(1) = rotationalVelocityMars*positionInertialFrame(0);
        velocityDueToMars(2) = 0.0;

        // Velocity from R-frame to I-frame
        const Eigen::Vector3d velocityInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(angleItoR)*velocityRotationalFrame+velocityDueToMars;

        /// Debug ///
//        std::cout<<"velocityRotationalFrame = "<<velocityRotationalFrame<<std::endl;
//        std::cout<<"velocityInertialFrame = "<<velocityInertialFrame<<std::endl;
        /// Debug ///

        tudat::basic_mathematics::Vector7d currentCartesianState;

        currentCartesianState(0) = positionInertialFrame(0); // x-position
        currentCartesianState(1) = positionInertialFrame(1); // y-position
        currentCartesianState(2) = positionInertialFrame(2); // z-position
        currentCartesianState(3) = velocityInertialFrame(0); // x-velocity
        currentCartesianState(4) = velocityInertialFrame(1); // y-velocity
        currentCartesianState(5) = velocityInertialFrame(2); // z-velocity
        currentCartesianState(6) = currentSphericalStateVector(6); // Mass


        outputVectorTSICartesian = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.

        // Filling the output vector
        outputVectorTSICartesian(0,0) = updatedStateAndTimeVector(7);   // Storing the updated time
        outputVectorTSICartesian(0,1) = currentCartesianState(0);   // Storing the updated x-position
        outputVectorTSICartesian(0,2) = currentCartesianState(1);   // Storing the updated y-position
        outputVectorTSICartesian(0,3) = currentCartesianState(2);   // Storing the updated z-position
        outputVectorTSICartesian(0,4) = currentCartesianState(3);   // Storing the updated x-velocity
        outputVectorTSICartesian(0,5) = currentCartesianState(4);   // Storing the updated y-velocity
        outputVectorTSICartesian(0,6) = currentCartesianState(5);   // Storing the updated z-velocity
        outputVectorTSICartesian(0,7) = currentCartesianState(6);   // Storing the updated MAV mass

        // Store the new cartesian values in the data storage matrix

        if (countTSI+countRKFTSI == 0){

          dataStoringMatrixTSICartesian.row(countTSI+countRKFTSI) = outputVectorTSICartesian.row(0); // Filling the matrix
        }
        else{
            dataStoringMatrixTSICartesian.conservativeResize(countTSI+countRKFTSI+1,8); // Making the matrix bigger in order to store more values

            dataStoringMatrixTSICartesian.row(countTSI+countRKFTSI) = outputVectorTSICartesian.row(0); // Filling the matrix
        }


//        std::cout<<"Current Cartesian State = "<<currentCartesianState<<std::endl;

        if (runningTimeTSI == 0.2){
                stateAtPoint2SecTSI = currentCartesianState;
                std::cout<<"The TSI state at 0.2 sec = "<<stateAtPoint2SecTSI<<std::endl;
        }

        currentStateAndTime.setCurrentStateAndTime(currentCartesianState,runningTimeTSI); // Update the current Cartesian state and time class!

        if (runningTimeTSI == coastStartTime){ // used to start the coasting phase
            MAV.setThrust(0); // Set the thrust equal to zero
        }

     countTSI++;

//     std::cout<<"countTSI = "<<countTSI<<std::endl;
//     }; // end of for-loop
//     std::cout<<"runningTimeTSI = "<<runningTimeTSI<<std::endl;
//     std::cout<<"endTimeTSI = "<<endTimeTSI<<std::endl;

    }while( !( endTimeTSI - runningTimeTSI <= std::numeric_limits< double >::epsilon( ) ) &&  !((currentSphericalStateAndTime.getCurrentSphericalRadius()-bodyReferenceRadius) >= EndAltitude));

        /// Adding the values to the file ///

        // Check if the file already exists.


        std::ifstream ifile2TSI(dataAbsolutePathTSI.c_str()); // Check it as an input file

        fexistsTSI = false;   // Set the default to "It does not exist"

        if (ifile2TSI){         // Attempt to open the file


           fexistsTSI = true;      // If the file can be opened it must exist

           ifile2TSI.close();   // Close the file

        }


        // If so: append, if not: error

        if (fexistsTSI == true){

            // Export the Taylor Series data matrix.
            std::ofstream exportFile1;                          // Define the file as an output file


            exportFile1.open(dataAbsolutePathTSI.c_str(),std::ios_base::app);      // Open the file in append mode

            exportFile1 << "\n";                                            // Make sure the new matrix start on a new line

            exportFile1 << dataStoringMatrixTSI.format( csvFormatTSI ); // Add the new values

            std::cout<<"The file called "<<dataAbsolutePathTSI<<" has been appended"<<std::endl;


            exportFile1.close( );   // Close the file

            // Export the Taylor Series "Spherical" Cartesian data matrix.
            std::ofstream exportFile1Cart;                          // Define the file as an output file


            exportFile1Cart.open(dataAbsolutePathTSICart.c_str(),std::ios_base::app);      // Open the file in append mode

            exportFile1Cart << "\n";                                            // Make sure the new matrix start on a new line

            exportFile1Cart << dataStoringMatrixTSICartesian.format( csvFormatTSI ); // Add the new values

            std::cout<<"The file called "<<dataAbsolutePathTSICart<<" has been appended"<<std::endl;


            exportFile1Cart.close( );   // Close the file
}
            else{

            std::cerr<<"Error: values could not be stored because storage file does not exist"<<std::endl;
        };

        std::cout<<"The TSI final state = "<<currentStateAndTime.getCurrentState()<<std::endl;
        std::cout<<"The final time ="<<currentStateAndTime.getCurrentTime()<<std::endl;
        std::cout<<"countTSI = "<<countTSI<<std::endl;

        /// Determine the CPU time taken ///


        // Determine the CPU time
        const double TSICPUTime = clock();
        // Determine the CPU time

        // Determine the elapsed CPU time

        const double elapsedTSICPUTime = TSICPUTime-initialCPUTime;

        std::cout<<"The elapsed TSI CPU time = "<<elapsedTSICPUTime/CLOCKS_PER_SEC<<" sec"<<std::endl;
                    std::cout<<"The elapsed TSI CPU time in clocks = "<<elapsedTSICPUTime<<std::endl;
//                    std::cout<<"CLOCKS_PER_SEC = "<<CLOCKS_PER_SEC<<std::endl;

//                    std::cout<<"sin(pi) = "<<sin(tudat::mathematical_constants::LONG_PI)<<std::endl;
//                    std::cout<<"cos(pi/2) = "<<cos(tudat::mathematical_constants::LONG_PI/2)<<std::endl;

        std::cout<<"////////////////////////////////////////////////////////////////// End of TSI //////////////////////////////////////////////////////////////////"<<std::endl;
//*/



        //// Reset the values ////
        MAV.resetThrust();
        MAV.resetReferenceArea();
//        std::cout<<"Thrust has been reset"<<std::endl;
//        std::cout<<"MAV.Thrust() = "<<MAV.Thrust()<<std::endl;
//        std::cout<<"runningTimeTSI = "<<runningTimeTSI<<std::endl;
        //// Reset the values ////



////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// Testing the RKF and other higher order integrators //////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

                    /// Setting the data collection file for RKF and inserting the first values ///

                        // Set directory where output files will be stored. THIS REQUIRES THE COMPLETE PATH IN ORDER TO WORK!!
                        const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/01.integrationResults/RKF/";


                        // Set output format for matrix output.
                        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

                        // Set absolute path to file containing the data.
                        std::string dataAbsolutePath = outputDirectory + "test4FullIntegrationRKF7(8)stateAndTime.csv";
                        std::string dataAbsolutePathSpherical = outputDirectory + "test4FullIntegrationRKF7(8)stateAndTime.csv";


                        // Create a row vector for the storing of the data
                        Eigen::MatrixXd outputVector = Eigen::MatrixXd::Zero(1,8); // Create a row vector for the storing of the data
                        Eigen::MatrixXd outputVectorSpherical = Eigen::MatrixXd::Zero(1,8); // Create a row vector for the storing of the data

                        // Getting the initial conditions for storage
//                        const tudat::basic_mathematics::Vector7d initialState = currentStateAndTime.getCurrentState();
//                        const tudat::basic_mathematics::Vector7d initialState = aState;
                        const tudat::basic_mathematics::Vector7d initialStateRKF = finalStateFirstSecond; // Get the initial state from the first RKF integration

                        std::cout<<"initialStateRKF = "<<initialStateRKF<<std::endl;

                        // Filling the output vector
//                        outputVector(0,0) = currentStateAndTime.getCurrentTime();   // Storing the initial time
                        outputVector(0,0) = finalTimeFirstSecond;                 // Setting the initial time
                        outputVector(0,1) = initialStateRKF(0);   // Storing the initial x position
                        outputVector(0,2) = initialStateRKF(1);   // Storing the initial y position
                        outputVector(0,3) = initialStateRKF(2);   // Storing the initial z position
                        outputVector(0,4) = initialStateRKF(3);   // Storing the initial x velocity
                        outputVector(0,5) = initialStateRKF(4);   // Storing the initial y velocity
                        outputVector(0,6) = initialStateRKF(5);   // Storing the initial z velocity
                        outputVector(0,7) = initialStateRKF(6);   // Storing the initial MAV mass

                        // Filling the output vector
                        outputVectorSpherical(0,0) = 0;                 // Setting the initial time
                        outputVectorSpherical(0,1) = sphericalState(0);   // Storing the initial radius
                        outputVectorSpherical(0,2) = sphericalState(1);   // Storing the initial latitude
                        outputVectorSpherical(0,3) = sphericalState(2);   // Storing the initial longitude
                        outputVectorSpherical(0,4) = sphericalState(3);   // Storing the initial velocity
                        outputVectorSpherical(0,5) = sphericalState(4);   // Storing the initial flight-path angle
                        outputVectorSpherical(0,6) = sphericalState(5);   // Storing the initial azimuth angle
                        outputVectorSpherical(0,7) = sphericalState(6);   // Storing the initial MAV mass




                        // Storing the data

                        std::ifstream ifile(dataAbsolutePath.c_str()); // Check it as an input file
//                        std::ifstream ifileSpherical(dataAbsolutePathSpherical.c_str()); // Check it as an input file


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
                            std::string newFileNameSpherical = "backupSphericalRKFFileAtDateAndTime_" + ComputerTimeString + ".csv";


                        std::cerr<<"The file name that you have chosen already exists, a new file with name "<<newFileName<<" will be created to store the data for now"<<std::endl;

                        // Set new absolute path to file containing the data.
                        dataAbsolutePath = outputDirectory + newFileName;
                        dataAbsolutePathSpherical = outputDirectory + newFileNameSpherical;


                        // Export the data.
                        std::ofstream exportFile1( dataAbsolutePath.c_str( ) ); // Make the new file
                        std::cout<<"New file called "<<dataAbsolutePath<<" has been created"<<std::endl;
                        exportFile1 << outputVector.format( csvFormat );          // Store the new values
                        exportFile1.close( );   // Close the file

                        std::ofstream exportFile1Spherical( dataAbsolutePathSpherical.c_str( ) ); // Make the new file
                        std::cout<<"New file called "<<dataAbsolutePathSpherical<<" has been created"<<std::endl;
                        exportFile1Spherical << outputVectorSpherical.format( csvFormat );          // Store the new values
                        exportFile1Spherical.close( );   // Close the file


                    }
                            else{

                            // Export the data.
                            std::ofstream exportFile1( dataAbsolutePath.c_str( ) ); // Make the new file
                            std::cout<<"New file called "<<dataAbsolutePath<<" has been created"<<std::endl;
                            exportFile1 << outputVector.format( csvFormat );          // Store the new values
                            exportFile1.close( );   // Close the file

                            std::ofstream exportFile1Spherical( dataAbsolutePathSpherical.c_str( ) ); // Make the new file
                            std::cout<<"New file called "<<dataAbsolutePathSpherical<<" has been created"<<std::endl;
                            exportFile1Spherical << outputVectorSpherical.format( csvFormat );          // Store the new values
                            exportFile1Spherical.close( );   // Close the file
                        };



                    /// Testing with the state derivative function class


                    // Full complete test

                    ascentStateDerivativeFunctionClass stateDerivativeFunctionClass(Mars,MAV);     // Initialize the class

                    // Set the initial values for the flight-path angle and heading angle
                    stateDerivativeFunctionClass.setFlightPathAngleAndHeadingAngle(FlightPathAngle,HeadingAngle);

                    // Just using boost::bind

                    boost::function< tudat::basic_mathematics::Vector7d( const double, const tudat::basic_mathematics::Vector7d ) > stateDerivativeFunction // Using boost::function to create the passing function
                            =boost::bind( &ascentStateDerivativeFunctionClass::ascentStateDerivativeFunction, &stateDerivativeFunctionClass, _1, _2);       // Then using boost::bind to bind the function as per class stateDerivativeFunctionClass which requires two inputs:
                                                                                                                                                            // _1 which is the current time and _2 which is the current state


                 ///// Testing the implementation in the integrator ///

                    // Initial conditions.
//                    const double initialTime = currentStateAndTime.getCurrentTime();                            // Time.
//                    const double initialTime = 0.0;                            // Time. set for verification
                    const double initialTime = finalTimeFirstSecond; // Set initial time to be equal to the final time of the first RKF integration
                //    Eigen::VectorXd initialState = currentStateAndTime.getCurrentState(); // State: start with zero velocity at the origin.

//                    const double endTime = setEndTime;     // Using the same time as defined for TSI
                    const double endTime = runningTimeTSI;     // Using the same time as defined for TSI

                    double coastStartTimeRKF; // used for the initial RKF loop

                    if (coastStartTime>=setEndTime){
                        coastStartTimeRKF = endTime;

                    }
                    else {
                        coastStartTimeRKF = coastStartTime;
                        coast = true;
                    }



                    // Step-size settings.
                    // The minimum and maximum step-size are set such that the input data is fully accepted by the
                    // integrator, to determine the steps to be taken.
                    const double zeroMinimumStepSize = std::numeric_limits< double >::epsilon( );
                    const double infiniteMaximumStepSize = std::numeric_limits< double >::infinity( );
//                    double stepSizeRKF = chosenStepSize;          // Using the same initial step-size as defined for TSI
                    double stepSizeRKF = finalStepSizeFirstSecond;          // Using the same initial step-size as defined for TSI


                    // Tolerances.
                    const double relativeTolerance = chosenLocalErrorTolerance;     // 1e-14 is used by TSI, original setting was 1e-15
                    const double absoluteTolerance = chosenLocalErrorTolerance;     // 1e-14 is used by TSI, original setting was 1e-15

                    // For RKF7(8) a step-size of 0.2 is only used if the tolerances are 1e-3.... and is accepted till 1e-8
                    // For RKF4(5) a step-size of 0.2 is only used if the tolerances are 1e-7.... and is accepted till 1e-10
                    // For DP8(7) a step-size of 0.2 is only used if the tolerances are 1e-7.... and is accepted till 1e-9

                    /// Set the initial state to be equal to the final time and state of the first RKF integration
//                    initialState = finalStateFirstSecond;

                    // Set name for integrator
                    string method;

////////////////////////// Choose your weapon! I mean integrator... ///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Simply uncomment the integrator you want and make sure the others are commented away ///////////////////////////////////////////////////////////

//                    /// RungeKutta4 numerical integrator.
//                    std::cout<<"You have chosen RK4 as your integration method"<<std::endl;
//                    tudat::numerical_integrators::RungeKutta4IntegratorXd RK4integrator(
//                        stateDerivativeFunction, initialTime, initialState );

//                    // Integrate to the specified end time.
//                    Eigen::VectorXd RK4endState = RK4integrator.integrateTo( endTime, stepSizeRKF );

//                    std::cout<<"Count = "<<endTime/stepSizeRKF<<std::endl;
//                    std::cout<<"The RK4 end state is "<<RK4endState<<std::endl;
//                    /// End of RungeKutta4 numerical integrator
///*                    // Is needed to comment the rest of the code since that is only used for the variable step-size methods



                    /// Runge-Kutta-Fehlberg 7(8) integrator.
                    std::cout<<"You have chosen RKF7(8) as your integration method"<<std::endl;
                    method = "RKF7(8)";
                       tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
                                   tudat::numerical_integrators::RungeKuttaCoefficients::get(
                                       tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78),
                                   stateDerivativeFunction, initialTime, initialStateRKF, zeroMinimumStepSize,
                                   infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

//                       /// Runge-Kutta-Fehlberg 4(5) integrator.
//                       std::cout<<"You have chosen RKF4(5) as your integration method"<<std::endl;
//                       method = "RKF4(5)";
//                          tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
//                                      tudat::numerical_integrators::RungeKuttaCoefficients::get(
//                                          tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg45),
//                                      stateDerivativeFunction, initialTime, initialStateRKF, zeroMinimumStepSize,
//                                      infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

//                          /// Dormand-Prince 8(7) integrator.
//                          std::cout<<"You have chosen DOPRIN8(7) as your integration method"<<std::endl;
//                          method = "DOPRIN8(7)";
//                             tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
//                                         tudat::numerical_integrators::RungeKuttaCoefficients::get(
//                                             tudat::numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince),
//                                         stateDerivativeFunction, initialTime, initialStateRKF, zeroMinimumStepSize,
//                                         infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                       std::cout<<"////////////////////////////////////////////////////////////////// Start of RKF //////////////////////////////////////////////////////////////////"<<std::endl;

                      // Set initial running time. This is updated after each step that the numerical integrator takes.
//                    double runningTime = 0.0;
                       double runningTime = initialTime; // Used to indicate that it has to start after the first RKF integration
                    int count = 0;


                    // Define storing matrix for the intermediate values
                    Eigen::MatrixXd dataStoringMatrix(1,8); // The size of this matrix will change in the do-loop
                    Eigen::MatrixXd dataStoringMatrixSpherical(1,8); // The size of this matrix will change in the do-loop
                    tudat::basic_mathematics::Vector7d stateAtPoint2SecRKF; // Storing the 0.2 seconds value specifically for comparison

                    bool RKFtimeCPUset = false; // Used to determine if the CPU timing cycle has begon
                    double initialCPUTimeRKF = 0.0; // Setting the default to zero

                    do
                    {
                        // Make sure the integrator does not integrate beyond the end time.

                        // Determine the CPU time
                        if (runningTime > 1.0 && RKFtimeCPUset == false){
                        initialCPUTimeRKF = clock();
                            RKFtimeCPUset = true;
}
                        // Determine the CPU time


                        if ( std::fabs( coastStartTimeRKF - runningTime )
                             <= std::fabs( stepSizeRKF ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) )
                        {
                            stepSizeRKF = coastStartTimeRKF - integrator.getCurrentIndependentVariable( );
                        }

//                        double prevStepSize = stepSizeRKF;

//                         std::cout<<"The current stepSize is "<<prevStepSize<<" s"<<std::endl;

                        // Perform a single integration step. Then update the step-size and running time.
                        integrator.performIntegrationStep( stepSizeRKF );
                        stepSizeRKF = integrator.getNextStepSize( );
                        runningTime = integrator.getCurrentIndependentVariable( );

                        Eigen::VectorXd currentState = integrator.getCurrentState();

//                        std::cout<<"The current stepSize is "<<prevStepSize<<" s"<<std::endl;
//                        std::cout<<"The current running time is "<<runningTime<<std::endl;
//                        if (runningTime <= 1.0 ){
////                            std::cout<<"CurrentState = "<<currentState<<std::endl;

//                        }



                        if (runningTime == 0.2){
                            std::cout<<"State at time 0.2 = "<<currentState<<std::endl;
                            stateAtPoint2SecRKF = currentState;
//                            std::cout<<"Latitude = "<<atan2(currentState(1),currentState(0))<<std::endl;
//                            std::cout<<"FlightPathAngle = "<<-asin((currentState(3)*(-cos(atan2(currentState(1),currentState(0))))+currentState(4)*(-sin(atan2(currentState(1),currentState(0)))))/
//                                                                   (sqrt(currentState(3)*currentState(3)+currentState(4)*currentState(4)+currentState(5)*currentState(5))))<<std::endl;

                        }

                        /// Get the spherical state from the initial few seconds of RKF

                        const double angleItoR_2 = rotationalVelocityMars*(inertialFrameTime+runningTime)-primeMeridianAngle;





                        // Position

                        Eigen::Vector3d CartesianPositionInertialFrame = Eigen::Vector3d::Zero(1,3); // Define the cartesian position inertial frame
                        CartesianPositionInertialFrame(0) = currentState(0); // x
                        CartesianPositionInertialFrame(1) = currentState(1); // y
                        CartesianPositionInertialFrame(2) = currentState(2); // z


                        const Eigen::Vector3d CartesianPositionRotationalFrame = tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix(angleItoR_2)*CartesianPositionInertialFrame;



                        Eigen::MatrixXd updatedStateAndTimeVector = Eigen::MatrixXd::Zero(1,8); // Creating a temporary vector to hold the converted state

                        updatedStateAndTimeVector(0) = sqrt(currentState(0)*currentState(0)+currentState(1)*currentState(1)+currentState(2)*currentState(2)); // Radius

                        updatedStateAndTimeVector(1) = asin(CartesianPositionRotationalFrame(2)/updatedStateAndTimeVector(0)) ; // latitude (delta)

                        updatedStateAndTimeVector(2) = atan2(CartesianPositionRotationalFrame(1),CartesianPositionRotationalFrame(0)); // longitude (tau)




                        /// Debug ///

                        /// Debug ///



                        // Velocity

                        Eigen::Vector3d CartesianVelocityInertialFrame = Eigen::Vector3d::Zero(1,3); // Define the cartesian velocity inertial frame
                        CartesianVelocityInertialFrame(0) = currentState(3); // Vx
                        CartesianVelocityInertialFrame(1) = currentState(4); // Vy
                        CartesianVelocityInertialFrame(2) = currentState(5); // Vz

                        Eigen::Vector3d rotatingPlanetCorrection = Eigen::Vector3d::Zero(3);

                        rotatingPlanetCorrection(0) = -rotationalVelocityMars*currentState(1);
                        rotatingPlanetCorrection(1) = rotationalVelocityMars*currentState(0);
                        rotatingPlanetCorrection(2) = 0.0;


                        const Eigen::Vector3d CartesianVelocityRotationalFrame = tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix(angleItoR_2)*(CartesianVelocityInertialFrame-rotatingPlanetCorrection); // Cartesian velocity in the rotational frame

                        const Eigen::Vector3d CartesianVelocityVerticalFrame = tudat::reference_frames::getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(updatedStateAndTimeVector(2),updatedStateAndTimeVector(1))*CartesianVelocityRotationalFrame; // Cartesian velocity in the local vertical frame


                        updatedStateAndTimeVector(3) = sqrt((currentState(3)+rotationalVelocityMars*currentState(1))*(currentState(3)+rotationalVelocityMars*currentState(1))+(currentState(4)-rotationalVelocityMars*currentState(0))*(currentState(4)-rotationalVelocityMars*currentState(0))+currentState(5)*currentState(5)); // Ground velocity scalar



                        updatedStateAndTimeVector(4) = -asin(CartesianVelocityVerticalFrame(2)/updatedStateAndTimeVector(3)); // flight-path angle

                //        std::cout<<"flight-path angle = "<<updatedStateAndTimeVectorRKFTSI(4)<<std::endl;


                        updatedStateAndTimeVector(5) = atan2(CartesianVelocityVerticalFrame(1),CartesianVelocityVerticalFrame(0)); // azimuth angle

                //        std::cout<<"azimuth angle = "<<updatedStateAndTimeVectorRKFTSI(5)<<std::endl;

                //        std::cout<<" "<<std::endl;


                        updatedStateAndTimeVector(6) = currentState(6); // MAV mass

                        /// Debug ///



                        /// Debug ///



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

                                outputVectorSpherical = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.


                                // Filling the spherical output vector
                                outputVectorSpherical(0,0) = integrator.getCurrentIndependentVariable();   // Storing the updated time
                                outputVectorSpherical(0,1) = updatedStateAndTimeVector(0);   // Storing the updated radius
                                outputVectorSpherical(0,2) = updatedStateAndTimeVector(1);   // Storing the updated latitude
                                outputVectorSpherical(0,3) = updatedStateAndTimeVector(2);   // Storing the updated longitude
                                outputVectorSpherical(0,4) = updatedStateAndTimeVector(3);   // Storing the updated velocity
                                outputVectorSpherical(0,5) = updatedStateAndTimeVector(4);   // Storing the updated flight-path angle
                                outputVectorSpherical(0,6) = updatedStateAndTimeVector(5);   // Storing the updated azimuth angle
                                outputVectorSpherical(0,7) = updatedStateAndTimeVector(6);   // Storing the updated MAV mass

                        // Store the new values in the data storage matrix

                        if (count == 0){

                          dataStoringMatrix.row(count) = outputVector.row(0); // Filling the matrix
                          dataStoringMatrixSpherical.row(count) = outputVectorSpherical.row(0); // Filling the matrix

                        }
                        else{
                            dataStoringMatrix.conservativeResize(count+1,8); // Making the matrix bigger in order to store more values

                            dataStoringMatrix.row(count) = outputVector.row(0); // Filling the matrix

                            dataStoringMatrixSpherical.conservativeResize(count+1,8); // Making the matrix bigger in order to store more values

                            dataStoringMatrixSpherical.row(count) = outputVectorSpherical.row(0); // Filling the matrix
                        }

//                        std::cout<<"outputVectorSphericalRKF = "<<outputVectorSpherical<<std::endl;


                        count++;


                    }while( !( coastStartTimeRKF - runningTime <= std::numeric_limits< double >::epsilon( ) ) &&  !((outputVectorSpherical(1)-bodyReferenceRadius) >= EndAltitude) );

 /// Coasting loop

 if (coast == true){
     MAV.setThrust(0); // Set the thrust zero for the coasting phase




     /// Testing with the state derivative function class


     // Full complete test

     ascentStateDerivativeFunctionClass stateDerivativeFunctionClassCoast(Mars,MAV);     // Initialize the class


     // Just using boost::bind

     boost::function< tudat::basic_mathematics::Vector7d( const double, const tudat::basic_mathematics::Vector7d ) > stateDerivativeFunctionCoast // Using boost::function to create the passing function
             =boost::bind( &ascentStateDerivativeFunctionClass::ascentStateDerivativeFunction, &stateDerivativeFunctionClassCoast, _1, _2);       // Then using boost::bind to bind the function as per class stateDerivativeFunctionClass which requires two inputs:
                                                                                                                                             // _1 which is the current time and _2 which is the current state





     tudat::basic_mathematics::Vector7d newInitialState;

     newInitialState(0) = outputVector(1); // x-position
     newInitialState(1) = outputVector(2); // y-position
     newInitialState(2) = outputVector(3); // z-position
     newInitialState(3) = outputVector(4); // x-velocity
     newInitialState(4) = outputVector(5); // y-velocity
     newInitialState(5) = outputVector(6); // z-velocity
     newInitialState(6) = outputVector(7); // MAV mass

     runningTime = outputVector(0); // current time



                    ////////////////////////// Choose your coasting weapon! I mean integrator... ///////////////////////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////// Simply uncomment the integrator you want and make sure the others are commented away ///////////////////////////////////////////////////////////

                    //                    /// RungeKutta4 numerical integrator.
                    //                    std::cout<<"You have chosen RK4 as your integration method"<<std::endl;
                    //                    tudat::numerical_integrators::RungeKutta4IntegratorXd RK4integrator(
                    //                        stateDerivativeFunction, initialTime, initialState );

                    //                    // Integrate to the specified end time.
                    //                    Eigen::VectorXd RK4endState = RK4integrator.integrateTo( endTime, stepSizeRKF );

                    //                    std::cout<<"Count = "<<endTime/stepSizeRKF<<std::endl;
                    //                    std::cout<<"The RK4 end state is "<<RK4endState<<std::endl;
                    //                    /// End of RungeKutta4 numerical integrator
                    ///*                    // Is needed to comment the rest of the code since that is only used for the variable step-size methods



                                        /// Runge-Kutta-Fehlberg 7(8) integrator.
                                        std::cout<<"You have chosen RKF7(8) as your integration method"<<std::endl;
                                        method = "RKF7(8)";
                                           tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integratorCoast(
                                                       tudat::numerical_integrators::RungeKuttaCoefficients::get(
                                                           tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78),
                                                       stateDerivativeFunctionCoast, runningTime, newInitialState, zeroMinimumStepSize,
                                                       infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

                    //                       /// Runge-Kutta-Fehlberg 4(5) integrator.
                    //                       std::cout<<"You have chosen RKF4(5) as your integration method"<<std::endl;
                    //                       method = "RKF4(5)";
                    //                          tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
                    //                                      tudat::numerical_integrators::RungeKuttaCoefficients::get(
                    //                                          tudat::numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg45),
                    //                                      stateDerivativeFunction, initialTime, initialState, zeroMinimumStepSize,
                    //                                      infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

                    //                          /// Dormand-Prince 8(7) integrator.
                    //                          std::cout<<"You have chosen DOPRIN8(7) as your integration method"<<std::endl;
                    //                          method = "DOPRIN8(7)";
                    //                             tudat::numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator(
                    //                                         tudat::numerical_integrators::RungeKuttaCoefficients::get(
                    //                                             tudat::numerical_integrators::RungeKuttaCoefficients::rungeKutta87DormandPrince),
                    //                                         stateDerivativeFunction, initialTime, initialState, zeroMinimumStepSize,
                    //                                         infiniteMaximumStepSize, relativeTolerance, absoluteTolerance );

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




                                        do
                                        {
                                            // Make sure the integrator does not integrate beyond the end time.

                                            // Determine the CPU time
                                            if (runningTime > 1.0 && RKFtimeCPUset == false){
                                            initialCPUTimeRKF = clock();
                                                RKFtimeCPUset = true;
                    }
                                            // Determine the CPU time


                                            if ( std::fabs( endTime - runningTime )
                                                 <= std::fabs( stepSizeRKF ) * ( 1.0 + std::numeric_limits< double >::epsilon( ) ) )
                                            {
                                                stepSizeRKF = endTime - integratorCoast.getCurrentIndependentVariable( );
                                            }

                    //                        double prevStepSize = stepSizeRKF;

                    //                         std::cout<<"The current stepSize is "<<prevStepSize<<" s"<<std::endl;

                                            // Perform a single integration step. Then update the step-size and running time.
                                            integratorCoast.performIntegrationStep( stepSizeRKF );
                                            stepSizeRKF = integratorCoast.getNextStepSize( );
                                            runningTime = integratorCoast.getCurrentIndependentVariable( );

                                            Eigen::VectorXd currentState = integratorCoast.getCurrentState();

//                                            std::cout<<"The current stepSize is "<<prevStepSize<<" s"<<std::endl;
//                                            std::cout<<"The current running time is "<<runningTime<<std::endl;
                    //                        if (runningTime <= 1.0 ){
                    ////                            std::cout<<"CurrentState = "<<currentState<<std::endl;

                    //                        }



                                            if (runningTime == 0.2){
                                                std::cout<<"State at time 0.2 = "<<currentState<<std::endl;
                                                stateAtPoint2SecRKF = currentState;
                    //                            std::cout<<"Latitude = "<<atan2(currentState(1),currentState(0))<<std::endl;
                    //                            std::cout<<"FlightPathAngle = "<<-asin((currentState(3)*(-cos(atan2(currentState(1),currentState(0))))+currentState(4)*(-sin(atan2(currentState(1),currentState(0)))))/
                    //                                                                   (sqrt(currentState(3)*currentState(3)+currentState(4)*currentState(4)+currentState(5)*currentState(5))))<<std::endl;

                                            }

                                            /// Get the spherical state from the initial few seconds of RKF

                                            const double angleItoR_2 = rotationalVelocityMars*(inertialFrameTime+runningTime)-primeMeridianAngle;





                                            // Position

                                            Eigen::Vector3d CartesianPositionInertialFrame = Eigen::Vector3d::Zero(1,3); // Define the cartesian position inertial frame
                                            CartesianPositionInertialFrame(0) = currentState(0); // x
                                            CartesianPositionInertialFrame(1) = currentState(1); // y
                                            CartesianPositionInertialFrame(2) = currentState(2); // z


                                            const Eigen::Vector3d CartesianPositionRotationalFrame = tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix(angleItoR_2)*CartesianPositionInertialFrame;



                                            Eigen::MatrixXd updatedStateAndTimeVector = Eigen::MatrixXd::Zero(1,8); // Creating a temporary vector to hold the converted state

                                            updatedStateAndTimeVector(0) = sqrt(currentState(0)*currentState(0)+currentState(1)*currentState(1)+currentState(2)*currentState(2)); // Radius

                                            updatedStateAndTimeVector(1) = asin(CartesianPositionRotationalFrame(2)/updatedStateAndTimeVector(0)) ; // latitude (delta)

                                            updatedStateAndTimeVector(2) = atan2(CartesianPositionRotationalFrame(1),CartesianPositionRotationalFrame(0)); // longitude (tau)




                                            /// Debug ///

                                            /// Debug ///



                                            // Velocity

                                            Eigen::Vector3d CartesianVelocityInertialFrame = Eigen::Vector3d::Zero(1,3); // Define the cartesian velocity inertial frame
                                            CartesianVelocityInertialFrame(0) = currentState(3); // Vx
                                            CartesianVelocityInertialFrame(1) = currentState(4); // Vy
                                            CartesianVelocityInertialFrame(2) = currentState(5); // Vz

                                            Eigen::Vector3d rotatingPlanetCorrection = Eigen::Vector3d::Zero(3);

                                            rotatingPlanetCorrection(0) = -rotationalVelocityMars*currentState(1);
                                            rotatingPlanetCorrection(1) = rotationalVelocityMars*currentState(0);
                                            rotatingPlanetCorrection(2) = 0.0;


                                            const Eigen::Vector3d CartesianVelocityRotationalFrame = tudat::reference_frames::getInertialToPlanetocentricFrameTransformationMatrix(angleItoR_2)*(CartesianVelocityInertialFrame-rotatingPlanetCorrection); // Cartesian velocity in the rotational frame

                                            const Eigen::Vector3d CartesianVelocityVerticalFrame = tudat::reference_frames::getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(updatedStateAndTimeVector(2),updatedStateAndTimeVector(1))*CartesianVelocityRotationalFrame; // Cartesian velocity in the local vertical frame


                                            updatedStateAndTimeVector(3) = sqrt((currentState(3)+rotationalVelocityMars*currentState(1))*(currentState(3)+rotationalVelocityMars*currentState(1))+(currentState(4)-rotationalVelocityMars*currentState(0))*(currentState(4)-rotationalVelocityMars*currentState(0))+currentState(5)*currentState(5)); // Ground velocity scalar



                                            updatedStateAndTimeVector(4) = -asin(CartesianVelocityVerticalFrame(2)/updatedStateAndTimeVector(3)); // flight-path angle

                                    //        std::cout<<"flight-path angle = "<<updatedStateAndTimeVectorRKFTSI(4)<<std::endl;


                                            updatedStateAndTimeVector(5) = atan2(CartesianVelocityVerticalFrame(1),CartesianVelocityVerticalFrame(0)); // azimuth angle

                                    //        std::cout<<"azimuth angle = "<<updatedStateAndTimeVectorRKFTSI(5)<<std::endl;

                                    //        std::cout<<" "<<std::endl;


                                            updatedStateAndTimeVector(6) = currentState(6); // MAV mass

                                            /// Debug ///



                                            /// Debug ///



                                            /// Storing the values ///



                                                    outputVector = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.


                                                    // Filling the output vector
                                                    outputVector(0,0) = integratorCoast.getCurrentIndependentVariable();   // Storing the updated time
                                                    outputVector(0,1) = currentState(0);   // Storing the updated x position
                                                    outputVector(0,2) = currentState(1);   // Storing the updated y position
                                                    outputVector(0,3) = currentState(2);   // Storing the updated z position
                                                    outputVector(0,4) = currentState(3);   // Storing the updated x velocity
                                                    outputVector(0,5) = currentState(4);   // Storing the updated y velocity
                                                    outputVector(0,6) = currentState(5);   // Storing the updated z velocity
                                                    outputVector(0,7) = currentState(6);   // Storing the updated MAV mass

                                                    outputVectorSpherical = Eigen::MatrixXd::Zero(1,8); // Setting the output vector to zero again to be sure.


                                                    // Filling the spherical output vector
                                                    outputVectorSpherical(0,0) = integratorCoast.getCurrentIndependentVariable();   // Storing the updated time
                                                    outputVectorSpherical(0,1) = updatedStateAndTimeVector(0);   // Storing the updated radius
                                                    outputVectorSpherical(0,2) = updatedStateAndTimeVector(1);   // Storing the updated latitude
                                                    outputVectorSpherical(0,3) = updatedStateAndTimeVector(2);   // Storing the updated longitude
                                                    outputVectorSpherical(0,4) = updatedStateAndTimeVector(3);   // Storing the updated velocity
                                                    outputVectorSpherical(0,5) = updatedStateAndTimeVector(4);   // Storing the updated flight-path angle
                                                    outputVectorSpherical(0,6) = updatedStateAndTimeVector(5);   // Storing the updated azimuth angle
                                                    outputVectorSpherical(0,7) = updatedStateAndTimeVector(6);   // Storing the updated MAV mass

                                            // Store the new values in the data storage matrix

                                            if (count == 0){

                                              dataStoringMatrix.row(count) = outputVector.row(0); // Filling the matrix
                                              dataStoringMatrixSpherical.row(count) = outputVectorSpherical.row(0); // Filling the matrix

                                            }
                                            else{
                                                dataStoringMatrix.conservativeResize(count+1,8); // Making the matrix bigger in order to store more values

                                                dataStoringMatrix.row(count) = outputVector.row(0); // Filling the matrix

                                                dataStoringMatrixSpherical.conservativeResize(count+1,8); // Making the matrix bigger in order to store more values

                                                dataStoringMatrixSpherical.row(count) = outputVectorSpherical.row(0); // Filling the matrix
                                            }

                    //                        if (currentState(6)<= EndMass){ // used to start the coasting phase
                    //                            MAV.setThrust(0); // Set the thrust equal to zero
                    //                        }

                    //std::cout<<"Current count = "<<count<<std::endl;
//                                            std::cout<<"(outputVectorSpherical(1)-bodyReferenceRadius) = "<<(outputVectorSpherical(1)-bodyReferenceRadius)<<std::endl;
                                            count++;


                                        }while( !( endTime - runningTime <= std::numeric_limits< double >::epsilon( ) ) &&  !((outputVectorSpherical(1)-bodyReferenceRadius) >= EndAltitude) );

 } // End of coasting phase

                    /// Storing the values to the file ///

//                    std::cout<<"It goes here"<<std::endl;


                    // Check if the file already exists.


                    std::ifstream ifile2(dataAbsolutePath.c_str()); // Check it as an input file
//                    std::ifstream ifile2Spherical(dataAbsolutePathSpherical.c_str()); // Check it as an input file

                    fexists = false;   // Set the default to "It does not exist"

                    if (ifile2){         // Attempt to open the file


                       fexists = true;      // If the file can be opened it must exist

                       ifile2.close();   // Close the file

                    }


                    // If so: append, if not: create new file and put data in

                    if (fexists == true){

                        // Export the data matrix.
                        std::ofstream exportFile1;                          // Define the file as an output file


                        exportFile1.open(dataAbsolutePath.c_str(),std::ios_base::app);      // Open the file in append mode

                        exportFile1 << "\n";                                            // Make sure the new matrix start on a new line

                        exportFile1 << dataStoringMatrix.format( csvFormat ); // Add the new values

                        std::cout<<"The file called "<<dataAbsolutePath<<" has been appended"<<std::endl;


                        exportFile1.close( );   // Close the file

                        std::ofstream exportFile1Spherical;                          // Define the file as an output file


                        exportFile1Spherical.open(dataAbsolutePathSpherical.c_str(),std::ios_base::app);      // Open the file in append mode

                        exportFile1Spherical << "\n";                                            // Make sure the new matrix start on a new line

                        exportFile1Spherical << dataStoringMatrixSpherical.format( csvFormat ); // Add the new values

                        std::cout<<"The file called "<<dataAbsolutePathSpherical<<" has been appended"<<std::endl;
            }
                        else{

                        std::cerr<<"Error: values could not be stored because storage file does not exist"<<std::endl;
                    };






                    // The result of the integration.
//                    Eigen::VectorXd endState = integrator.getCurrentState( );
                    tudat::basic_mathematics::Vector7d endState;

                    endState(0) = outputVector(1); // x-position
                    endState(1) = outputVector(2); // y-position
                    endState(2) = outputVector(3); // z-position
                    endState(3) = outputVector(4); // x-velocity
                    endState(4) = outputVector(5); // y-velocity
                    endState(5) = outputVector(6); // z-velocity
                    endState(6) = outputVector(7); // MAV mass


//                    std::cout<<"Final number of integration steps is "<<count-countRKFTSI<<std::endl;
                    std::cout<<"Final number of integration steps is "<<count<<std::endl;
                    std::cout<<"The end state is "<<endState<<std::endl;

                    /// Determine the CPU time taken ///


                    // Determine the CPU time
                    const double finalCPUTime = clock();
                    // Determine the CPU time

                    // Determine the elapsed CPU time

                    const double elapsedCPUTime = finalCPUTime-initialCPUTimeRKF;

                    std::cout<<"The elapsed RKF CPU time = "<<elapsedCPUTime/CLOCKS_PER_SEC<<" sec"<<std::endl;
                    std::cout<<"The elapsed RKF CPU time in clocks = "<<elapsedCPUTime<<std::endl;
//                    std::cout<<"CLOCKS_PER_SEC = "<<CLOCKS_PER_SEC<<std::endl;

//                    std::cout<<"sin(pi) = "<<sin(tudat::mathematical_constants::LONG_PI)<<std::endl;
//                    std::cout<<"cos(pi/2) = "<<cos(tudat::mathematical_constants::LONG_PI/2)<<std::endl;

                    std::cout<<"////////////////////////////////////////////////////////////////// End of RKF //////////////////////////////////////////////////////////////////"<<std::endl;

                    //*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




                  /// Compute the differences between the methods ///

                    if (comparison == true){
                    // Difference t = 0.2 sec

//                    tudat::basic_mathematics::Vector7d differenceFractionPoint2Sec;
//                    tudat::basic_mathematics::Vector7d differencePoint2Sec;

////                    std::cout<<"stateAtPoint2SecTSI = "<<stateAtPoint2SecTSI<<std::endl;
////                    std::cout<<"stateAtPoint2SecRKF = "<<stateAtPoint2SecRKF<<std::endl;

//                    for (int i = 0;i<8;i++){
//                    differenceFractionPoint2Sec(i) = stateAtPoint2SecTSI(i)/stateAtPoint2SecRKF(i);
//                    differencePoint2Sec(i) = stateAtPoint2SecTSI(i)-stateAtPoint2SecRKF(i);
//}
//                    std::cout<<"The difference fraction between TSI and "<<method<<" for t = 0.2 sec = "<<"\n"<<
//                               differenceFractionPoint2Sec<<std::endl;
//                    std::cout<<"The difference between TSI and "<<method<<" for t = 0.2 sec = "<<"\n"<<
//                               differencePoint2Sec<<std::endl;

                    // Difference end state

                    tudat::basic_mathematics::Vector7d differenceFractionEnd;
                    tudat::basic_mathematics::Vector7d differenceEnd;
                    tudat::basic_mathematics::Vector7d endStateTSI = currentStateAndTime.getCurrentState();

//                    std::cout<<"endStateTSI = "<<endStateTSI<<std::endl;
//                    std::cout<<"endStateRKF = "<<endState<<std::endl;

                    for (int i = 0;i<8;i++){
                    differenceFractionEnd(i) = endStateTSI(i)/endState(i);
                    differenceEnd(i) = endStateTSI(i)-endState(i);
}
                    std::cout<<"The difference fraction between TSI and "<<method<<" at the end time = "<<"\n"<<
                               differenceFractionEnd<<std::endl;
                    std::cout<<"The difference between TSI and "<<method<<" at the end time = "<<"\n"<<
                               differenceEnd<<std::endl;

                    std::cout<<"/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
                    } // end comparison



//*/

///////////////////////////////////////////////////////////////////// Circularisation ///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /// Convert the inertial Cartesian coordinates for both methods to Kepler elements ///
//    keplerianElements( 0 ) = semiMajorAxis,                                             [km]
//    keplerianElements( 1 ) = eccentricity,                                              [-]
//    keplerianElements( 2 ) = inclination,                                             [rad]
//    keplerianElements( 3 ) = argument of periapsis,                                   [rad]
//    keplerianElements( 4 ) = longitude of ascending node,                             [rad]
//    keplerianElements( 5 ) = true anomaly.                                            [rad]



       tudat::basic_mathematics::Vector6d RKFendCartesianCoordinates; // The final conditions from the RKF integrator

       RKFendCartesianCoordinates(0) = endState(0); // x-position
       RKFendCartesianCoordinates(1) = endState(1); // y-position
       RKFendCartesianCoordinates(2) = endState(2); // z-position
       RKFendCartesianCoordinates(3) = endState(3); // x-velocity
       RKFendCartesianCoordinates(4) = endState(4); // y-velocity
       RKFendCartesianCoordinates(5) = endState(5); // z-velocity

       tudat::basic_mathematics::Vector7d TSIendState = currentStateAndTime.getCurrentState(); // This TSI integrator end state includes the MAV mass (which is not required for the Kepler Elements)

       tudat::basic_mathematics::Vector6d TSIendCartesianCoordinates;  // The final conditions from the TSI integrator

       TSIendCartesianCoordinates(0) = TSIendState(0); // x-position
       TSIendCartesianCoordinates(1) = TSIendState(1); // y-position
       TSIendCartesianCoordinates(2) = TSIendState(2); // z-position
       TSIendCartesianCoordinates(3) = TSIendState(3); // x-velocity
       TSIendCartesianCoordinates(4) = TSIendState(4); // y-velocity
       TSIendCartesianCoordinates(5) = TSIendState(5); // z-velocity


//         Convert the integrator end coordinates
        tudat::basic_mathematics::Vector6d  RKFendKeplerElements = tudat::orbital_element_conversions::convertCartesianToKeplerianElements(RKFendCartesianCoordinates,Mars.standardGravitationalParameter()); // RKF
        tudat::basic_mathematics::Vector6d  TSIendKeplerElements = tudat::orbital_element_conversions::convertCartesianToKeplerianElements(TSIendCartesianCoordinates,Mars.standardGravitationalParameter()); // TSI

//        std::cout<<"RKFendKeplerElements = "<<"\n"<<RKFendKeplerElements<<std::endl;
//        std::cout<<"TSIendKeplerElements = "<<"\n"<<TSIendKeplerElements<<std::endl;

        // Orbital velocities

        const double currentRKForbitalVelocity = sqrt(Mars.standardGravitationalParameter()*(2/(sqrt(endState(0)*endState(0)+endState(1)*endState(1)+endState(2)*endState(2)))-1/RKFendKeplerElements(0))); // V_RKF
        const double currentTSIorbitalVelocity = sqrt(Mars.standardGravitationalParameter()*(2/(sqrt(TSIendState(0)*TSIendState(0)+TSIendState(1)*TSIendState(1)+TSIendState(2)*TSIendState(2)))-1/TSIendKeplerElements(0))); // V_TSI


//        std::cout<<"currentRKForbitalVelocity = "<<currentRKForbitalVelocity<<std::endl;
//        std::cout<<"currentTSIorbitalVelocity = "<<currentTSIorbitalVelocity<<std::endl;

        const double requiredOrbitalVelocity = sqrt(Mars.standardGravitationalParameter()/(bodyReferenceRadius+EndAltitude));

//        std::cout<<"requiredOrbitalVelocity = "<<requiredOrbitalVelocity<<std::endl;

//        std::cout<<"currentRKFradius = "<<sqrt(endState(0)*endState(0)+endState(1)*endState(1)+endState(2)*endState(2))<<std::endl;
//        std::cout<<"currentTSIradius = "<<sqrt(TSIendState(0)*TSIendState(0)+TSIendState(1)*TSIendState(1)+TSIendState(2)*TSIendState(2))<<std::endl;
//        std::cout<<"requiredRadius = "<<(bodyReferenceRadius+EndAltitude)<<std::endl;
//        std::cout<<"currentRKFpericentreRadius = "<<RKFendKeplerElements(0)*(1.0-RKFendKeplerElements(1))<<std::endl;
//        std::cout<<"currentTSIpericentreRadius = "<<TSIendKeplerElements(0)*(1.0-TSIendKeplerElements(1))<<std::endl;
//        std::cout<<"currentRKFapocentreRadius = "<<RKFendKeplerElements(0)*(1.0+RKFendKeplerElements(1))<<std::endl;
//        std::cout<<"currentTSIapocentreRadius = "<<TSIendKeplerElements(0)*(1.0+TSIendKeplerElements(1))<<std::endl;

//        std::cout<<"currentInertialVelocity = "<<sqrt(endState(3)*endState(3)+endState(4)*endState(4)+endState(5)*endState(5))<<std::endl;




        // Flight path angles

        const double pRKF = RKFendKeplerElements(0)*(1.0-RKFendKeplerElements(1)*RKFendKeplerElements(1)); // Semi-latus rectum for RKF
        const double cosRKFflightPathAngle = (pRKF/sqrt(endState(0)*endState(0)+endState(1)*endState(1)+endState(2)*endState(2)))/
                (sqrt(2*(pRKF/sqrt(endState(0)*endState(0)+endState(1)*endState(1)+endState(2)*endState(2)))-
                      (1.0-RKFendKeplerElements(1)*RKFendKeplerElements(1)))); // Cosine of the orbital flight-path angle for RKF
//        const double RKFflightPathAngle = acos(cosRKFflightPathAngle); // Orbital flight-path angle for RKF


        const double pTSI = TSIendKeplerElements(0)*(1.0-TSIendKeplerElements(1)*TSIendKeplerElements(1)); // Semi-latus rectum for TSI
        const double cosTSIflightPathAngle = (pTSI/sqrt(TSIendState(0)*TSIendState(0)+TSIendState(1)*TSIendState(1)+TSIendState(2)*TSIendState(2)))/
                (sqrt(2*(pTSI/sqrt(TSIendState(0)*TSIendState(0)+TSIendState(1)*TSIendState(1)+TSIendState(2)*TSIendState(2)))-
                      (1.0-TSIendKeplerElements(1)*TSIendKeplerElements(1)))); // Cosine of the orbital flight-path angle for TSI
//        const double TSIflightPathAngle = acos(cosTSIflightPathAngle); // Orbital flight-path angle for TSI

//        std::cout<<"RKFflightPathAngle = "<<RKFflightPathAngle<<std::endl;
//        std::cout<<"TSIflightPathAngle = "<<TSIflightPathAngle<<std::endl;

//        std::cout<<"testRKFVelocity = "<<sqrt(requiredOrbitalVelocity*requiredOrbitalVelocity*(2.0-(sqrt(endState(0)*endState(0)+endState(1)*endState(1)+endState(2)*endState(2))/pRKF)*
//                                                                                               (1.0-RKFendKeplerElements(1)*RKFendKeplerElements(1))))<<std::endl;

        // Delta-V's

//        const double deltaVforRKF = sqrt(requiredOrbitalVelocity*requiredOrbitalVelocity+currentRKForbitalVelocity*currentRKForbitalVelocity-
//                                         2*requiredOrbitalVelocity*currentRKForbitalVelocity*cosRKFflightPathAngle); // Delta-V required for the RKF change into the circular orbit

        const double deltaVforRKF = sqrt(requiredOrbitalVelocity*requiredOrbitalVelocity+currentRKForbitalVelocity*currentRKForbitalVelocity-
                                         2*requiredOrbitalVelocity*currentRKForbitalVelocity*cosRKFflightPathAngle*cos(desiredInclination-RKFendKeplerElements(2))); // Delta-V required for the RKF change into the circular orbit including inclination change

//        const double deltaVforTSI = sqrt(requiredOrbitalVelocity*requiredOrbitalVelocity+currentTSIorbitalVelocity*currentTSIorbitalVelocity-
//                                         2*requiredOrbitalVelocity*currentTSIorbitalVelocity*cosTSIflightPathAngle); // Delta-V required for the TSI change into the circular orbit

        const double deltaVforTSI = sqrt(requiredOrbitalVelocity*requiredOrbitalVelocity+currentTSIorbitalVelocity*currentTSIorbitalVelocity-
                                         2*requiredOrbitalVelocity*currentTSIorbitalVelocity*cosTSIflightPathAngle*cos(desiredInclination-TSIendKeplerElements(2))); // Delta-V required for the TSI change into the circular orbit including inclination change

        std::cout<<"deltaVforRKF = "<<deltaVforRKF<<std::endl;
        std::cout<<"deltaVforTSI = "<<deltaVforTSI<<std::endl;

        std::cout<<"desiredInclination-RKFendKeplerElements(2) = "<<desiredInclination-RKFendKeplerElements(2)<<std::endl;
        std::cout<<"desiredInclination-TSIendKeplerElements(2) = "<<desiredInclination-TSIendKeplerElements(2)<<std::endl;

        std::cout<<"desiredInclination = "<<rad2deg(desiredInclination)<<std::endl;
        std::cout<<"RKFInclination = "<<rad2deg(RKFendKeplerElements(2))<<std::endl;
        std::cout<<"TSIInclination = "<<rad2deg(TSIendKeplerElements(2))<<std::endl;




        // Required propellant mass

        const double propMassRKF = endState(6)-endState(6)/(exp(deltaVforRKF/(MAV.specificImpulse()*tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION))); // RKF propellant mass required (minimum)
        const double propMassTSI = TSIendState(6)-TSIendState(6)/(exp(deltaVforTSI/(MAV.specificImpulse()*tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION))); // TSI propellant mass required (minimum)


        std::cout<<"propMassRKF = "<<propMassRKF<<std::endl;
        std::cout<<"propMassTSI = "<<propMassTSI<<std::endl;


        // Final MAV mass

        const double finalMassRKF = endState(6)-propMassRKF; // RKF final MAV mass
        const double finalMassTSI = TSIendState(6)-propMassTSI; // TSI final MAV mass

        std::cout<<"finalMassRKF = "<<finalMassRKF<<std::endl;
        std::cout<<"finalMassTSI = "<<finalMassTSI<<std::endl;



    return 0;
}


