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
#include <iomanip>
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

// testing


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


int main()

{

std::cout<<setprecision(15)<<"Setting output precision to 15"<<std::endl;

    /// Testing the Celestial Body class ///


//    // First test

//    celestialBody Mars;

    // Second test

//    const std::string planet = "Mars";
//    const std::string planet = "Venus";

    celestialBody Mars;

//    celestialBody Mars(planet);
//    Mars.setPlanet(planet);  // Does not exist in the class anymore!

   const double adiabeticIndex = Mars.adiabeticIndex();
   const double specificGasConstant = Mars.specificGasConstant();
   const double standardGravitationalParameter = Mars.standardGravitationalParameter();
   const double rotationalVelocity = Mars.rotationalVelocity();
   const double primeMeridianAngle = Mars.primeMeridianAngle();
   const double inertialFrameTime = Mars.inertialFrameTime();

  const double bodyReferenceRadius = Mars.bodyReferenceRadius();


   const Eigen::MatrixXd temperaturePolyCoefficients = Mars.temperaturePolyCoefficients();
   const Eigen::MatrixXd temperatureAltitudeRanges = Mars.temperatureAltitudeRanges();
   const Eigen::VectorXd densityPolyCoefficients = Mars.densityPolyCoefficients();
/*
       std::cout<<"The adiabeticIndex for the Martian atmosphere is "<<adiabeticIndex<<std::endl;
       std::cout<<"The specificGasConstant for the Martian atmosphere is "<<specificGasConstant<<std::endl;
       std::cout<<"The standardGravitationalParameter for the Martian atmosphere is "<<standardGravitationalParameter<<std::endl;
       std::cout<<"The rotationalVelocity for the Martian atmosphere is "<<rotationalVelocity<<std::endl;
       std::cout<<"The primeMeridianAngle for the Martian atmosphere is "<<primeMeridianAngle<<std::endl;
       std::cout<<"The inertialFrameTime for the Martian atmosphere is "<<inertialFrameTime<<std::endl;
       std::cout<<"The marsMolaRadius for the Martian atmosphere is "<<bodyReferenceRadius<<std::endl;
       std::cout<<"The temperaturePolyCoefficients for the Martian atmosphere is "<<temperaturePolyCoefficients<<std::endl;
       std::cout<<"The temperatureAltitudeRanges for the Martian atmosphere is "<<temperatureAltitudeRanges<<std::endl;
       std::cout<<"The densityPolyCoefficients for the Martian atmosphere is "<<densityPolyCoefficients<<std::endl;



//   boost::shared_ptr< celestialBody > MarsTwo = boost::make_shared< celestialBody > ();     // Still not sure how this works so lets keep it simple for now shall we...

//     boost::shared_ptr< Mars.adiabeticIndex() > adiabeticIndexTwoTest = boost::make_shared< Mars.adiabeticIndex() > ();

//            std::cout<<"The adiabeticIndexTwoTest for the Martian atmosphere is "<<adiabeticIndexTwoTest<<std::endl;

*/
    /// Testing the vehicle class ///

    MarsAscentVehicle MAV;

    // Returning the different constant parameters

    const double Thrust = MAV.Thrust();                                                          // T     engine nominal thrust
    const double specificImpulse = MAV.specificImpulse();                                        // Isp     engine nominal specific impulse
    const double referenceArea = MAV.referenceArea();                                             // S   vehicle reference area
    const Eigen::MatrixXd dragCoefficientPolyCoefficients = MAV.dragCoefficientPolyCoefficients();            // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
    const Eigen::MatrixXd dragCoefficientMachRanges = MAV.dragCoefficientMachRanges();                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient
    const Eigen::MatrixXd thrustAzimuth = MAV.thrustAzimuth();                                // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
    const Eigen::MatrixXd thrustElevation = MAV.thrustElevation();                                // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)

/*
    std::cout<<"The Thrust for the Mars Ascent Vehicle is "<<Thrust<<std::endl;
    std::cout<<"The specificImpulse for the Mars Ascent Vehicle is "<<specificImpulse<<std::endl;
    std::cout<<"The referenceArea for the Mars Ascent Vehicle is "<<referenceArea<<std::endl;
    std::cout<<"The dragCoefficientPolyCoefficients for the Mars Ascent Vehicle is "<<dragCoefficientPolyCoefficients<<std::endl;
    std::cout<<"The dragCoefficientMachRanges for the Mars Ascent Vehicle is "<<dragCoefficientMachRanges<<std::endl;
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

/*
     Initial conditions

    tudat::basic_mathematics::Vector6d aState;

    aState(0) = 1;
    aState(1) = 1;
    aState(2) = 1;
    aState(3) = 2;
    aState(4) = 2;
    aState(5) = 2;

    double aMass = 227;   // [kg] from literature study


    /// Test original Body class converted to StateAndTime class /// (The pointer doesn't really do much, is too much work and end up with the same thing in the end. Only useful if it actually calls a function with input variables.)


//    const satellite_propagator_examples::StateAndTimePointer currentStateAndTime = boost::make_shared< satellite_propagator_examples::StateAndTime > (aState,aMass); // Creating the current state class using the namespace by pointing

//    boost::function<const tudat::basic_mathematics::Vector6d ()> currentState_ = boost::bind( &satellite_propagator_examples::StateAndTime::getCurrentState, currentStateAndTime);
//    boost::function<const Eigen::Vector3d ()> currentPosition_ = boost::bind( &satellite_propagator_examples::StateAndTime::getCurrentPosition, currentStateAndTime);
//    boost::function<const Eigen::Vector3d ()> currentVelocity_ = boost::bind( &satellite_propagator_examples::StateAndTime::getCurrentVelocity, currentStateAndTime);
//    boost::function<const double ()> currentMass_ = boost::bind( &satellite_propagator_examples::StateAndTime::getCurrentMass, currentStateAndTime);
//    boost::function<const double ()> currentTime_ = boost::bind( &satellite_propagator_examples::StateAndTime::getCurrentTime, currentStateAndTime);


//        const tudat::basic_mathematics::Vector6d currentState = currentState_();
//        const Eigen::Vector3d currentPosition = currentPosition_();
//        const Eigen::Vector3d currentVelocity = currentVelocity_();
//        const double currentMass = currentMass_();
//        const double currentTime = currentTime_();


//    satellite_propagator_examples::StateAndTime currentStateAndTime(aState,aMass);        // Creating the current state class using the namespace and class directly

//    const tudat::basic_mathematics::Vector6d currentState = currentStateAndTime.getCurrentState();
//    const Eigen::Vector3d currentPosition = currentStateAndTime.getCurrentPosition();
//    const Eigen::Vector3d currentVelocity = currentStateAndTime.getCurrentVelocity();
//    const double currentMass = currentStateAndTime.getCurrentMass();
//    const double currentTime = currentStateAndTime.getCurrentTime();



    /// Test modified StateAndTime class (removed unnecessary namespace) /// (Still, pointer is useless I think)



//        const StateAndTimePointer currentStateAndTime = boost::make_shared< StateAndTime > (aState,aMass); // Creating the current state class using the namespace by pointing

//        boost::function<const tudat::basic_mathematics::Vector6d ()> currentState_ = boost::bind( &StateAndTime::getCurrentState, currentStateAndTime);
//        boost::function<const Eigen::Vector3d ()> currentPosition_ = boost::bind( &StateAndTime::getCurrentPosition, currentStateAndTime);
//        boost::function<const Eigen::Vector3d ()> currentVelocity_ = boost::bind( &StateAndTime::getCurrentVelocity, currentStateAndTime);
//        boost::function<const double ()> currentMass_ = boost::bind( &StateAndTime::getCurrentMass, currentStateAndTime);
//        boost::function<const double ()> currentTime_ = boost::bind( &StateAndTime::getCurrentTime, currentStateAndTime);


//            const tudat::basic_mathematics::Vector6d currentState = currentState_();
//            const Eigen::Vector3d currentPosition = currentPosition_();
//            const Eigen::Vector3d currentVelocity = currentVelocity_();
//            const double currentMass = currentMass_();
//            const double currentTime = currentTime_();




//        StateAndTime currentStateAndTime(aState,aMass);        // Creating the current state class using the namespace and class directly

//        const tudat::basic_mathematics::Vector6d currentState = currentStateAndTime.getCurrentState();
//        const Eigen::Vector3d currentPosition = currentStateAndTime.getCurrentPosition();
//        const Eigen::Vector3d currentVelocity = currentStateAndTime.getCurrentVelocity();
//        const double currentMass = currentStateAndTime.getCurrentMass();
//        const double currentTime = currentStateAndTime.getCurrentTime();

*/



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


//    std::cout<<setprecision(15)<<"initialLatitude = "<<initialLatitude<<std::endl;

    const double initialRadius = bodyReferenceRadius+initialAltitude;               // Starting radius in km
//        const double initialRadius = bodyReferenceRadius+initialAltitude;               // Starting radius in m
/*
    std::cout<<"initialRadius = "<<initialRadius<<std::endl;
    std::cout<<"initialAltitude = "<<initialAltitude<<std::endl;
*/
/*
//    const Eigen::Vector3d initialSphericalPosition = Eigen::Vector3d(initialRadius,initialLatitude,initialLongitude);

//        std::cout<<"The initialSpehricalPosition is "<<initialSphericalPosition<<std::endl;
*/
        // Converting the initial spherical position to cartesian position using the standard convertSphericalToCartesian function of Tudat
        // Please note that this function requires the zenith angle as input which is pi/2-latitude!

//    const Eigen::Vector3d initialCartesianPositionRotationalFrame = tudat::coordinate_conversions::convertSphericalToCartesian(Eigen::Vector3d(initialRadius,((tudat::mathematical_constants::LONG_PI/2)-initialLatitude),initialLongitude));

      Eigen::Vector3d initialCartesianPositionRotationalFrame = Eigen::Vector3d::Zero(3);

      initialCartesianPositionRotationalFrame(0) = initialRadius*cos(initialLatitude)*cos(initialLongitude); // x_R
      initialCartesianPositionRotationalFrame(1) = initialRadius*cos(initialLatitude)*sin(initialLongitude); // y_R
      initialCartesianPositionRotationalFrame(2) = initialRadius*sin(initialLatitude); // z_R


//    const Eigen::Vector3d initialCartesianPositionRotationalFrame = tudat::coordinate_conversions::convertSphericalToCartesian(Eigen::Vector3d(initialRadius,((3.1415926535897932384626433832795028841971693993751058/2)-initialLatitude),initialLongitude)); // This also made it crash..


//    const Eigen::Vector3d initialCartesianPositionRotationalFrame = tudat::coordinate_conversions::convertSphericalToCartesian(Eigen::Vector3d(initialRadius,((tudat::mathematical_constants::PI/2)-initialLatitude),initialLongitude)); // This makes it crash and say: The program has unexpectedly finished.

//    const Eigen::Vector3d initialCartesianPositionInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocity*inertialFrameTime-primeMeridianAngle)*initialCartesianPositionRotationalFrame;


    const Eigen::Vector3d initialCartesianPositionInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocity*inertialFrameTime-primeMeridianAngle)*initialCartesianPositionRotationalFrame;




    // Compute initial velocity in y-direction as seen from the launch site in the inertial frame

    const Eigen::Vector3d initialVelocityLaunchSite = Eigen::Vector3d(0,(rotationalVelocity*initialRadius*cos(initialLatitude)),0);

//    const Eigen::Vector3d initialVelocityInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocity*inertialFrameTime-primeMeridianAngle+initialLongitude)*initialVelocityLaunchSite;


    const Eigen::Vector3d initialVelocityInertialFrame = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(rotationalVelocity*inertialFrameTime-primeMeridianAngle+initialLongitude)*initialVelocityLaunchSite;


/*    /// Debug ///

    std::cout<<"x1 main = "<<initialCartesianPositionInertialFrame(0)<<std::endl;
    std::cout<<"x2 main = "<<initialCartesianPositionInertialFrame(1)<<std::endl;
    std::cout<<"x1 main-847113.311014168 = "<<initialCartesianPositionInertialFrame(0)-847113.311014168<<std::endl;
    std::cout<<"x2 main-3054591.91823781 = "<<initialCartesianPositionInertialFrame(1)-3054591.91823781<<std::endl;
*/


/*
        std::cout<<initialLatitude<<std::endl;
        std::cout<<initialLongitude<<std::endl;
        std::cout<<"The initialCartesianPositionRotationalFrame is "<<initialCartesianPositionRotationalFrame<<std::endl;
        std::cout<<"The initialCartesianPositionInertialFrame is "<<initialCartesianPositionInertialFrame<<std::endl;
        std::cout<<"The initialVelocityLaunchSite is "<<initialVelocityLaunchSite<<std::endl;
        std::cout<<"The initialVelocityInertialFrame is "<<initialVelocityInertialFrame<<std::endl;

*/
        /// Testing StateAndTime class using modified vector ///

    tudat::basic_mathematics::Vector7d aState;

    aState(0) = initialCartesianPositionInertialFrame(0);
    aState(1) = initialCartesianPositionInertialFrame(1);
    aState(2) = initialCartesianPositionInertialFrame(2);
    aState(3) = initialVelocityInertialFrame(0);
    aState(4) = initialVelocityInertialFrame(1);
    aState(5) = initialVelocityInertialFrame(2);
    aState(6) = 227;  // Mass [kg] from literature study
/*
    std::cout<<"initialRadius minus computed initial position inertial = "<<initialRadius-sqrt(aState(0)*aState(0)+aState(1)*aState(1)+aState(2)*aState(2))<<std::endl;
    std::cout<<"initialRadius minus computed initial position rotational = "<<initialRadius-sqrt(initialCartesianPositionRotationalFrame(0)*initialCartesianPositionRotationalFrame(0)+
                                                                                                 initialCartesianPositionRotationalFrame(1)*initialCartesianPositionRotationalFrame(1)+
                                                                                                 initialCartesianPositionRotationalFrame(2)*initialCartesianPositionRotationalFrame(2))<<std::endl;
*/

    StateAndTime currentStateAndTime(aState);        // Creating the current state class using the namespace and class directly

//    const tudat::basic_mathematics::Vector7d currentState = currentStateAndTime.getCurrentState();
//    const Eigen::Vector3d currentPosition = currentStateAndTime.getCurrentPosition();
//    const Eigen::Vector3d currentVelocity = currentStateAndTime.getCurrentVelocity();
    const double currentMass = currentStateAndTime.getCurrentMass();
    const double currentTime = currentStateAndTime.getCurrentTime();


/*
    std::cout<<"The currentState is "<<currentState<<std::endl;
    std::cout<<"The currentPosition is "<<currentPosition<<std::endl;
    std::cout<<"The currentVelocity is "<<currentVelocity<<std::endl;
    std::cout<<"The currentMass is "<<currentMass<<std::endl;
    std::cout<<"The currentTime is "<<currentTime<<std::endl;
*/



    /// Thrust acceleration in B-frame ///   thrustAccelerationsBframe

    const Eigen::Vector3d thrustAccelerationsPframe = Eigen::Vector3d((Thrust/currentMass),0,0);

    const double thrustAzimuthTestDeg = 0;             // thrust azimuth gimbal angle [Deg] 10 for testing
    const double thrustElevationTestDeg = 0;            // thrust elevation gimbal angle [Deg] 5 for testing

    const double thrustAzimuthTest = deg2rad(thrustAzimuthTestDeg);     // thrust azimuth gimbal angle [rad]
    const double thrustElevationTest = deg2rad(thrustElevationTestDeg); // thrust elevation gimbal angle [rad]

//    const Eigen::Vector3d thrustAccelerationsBframe = tudat::reference_frames::getPropulsionToBodyFrameTransformationMatrix(thrustAzimuthTest,thrustElevationTest)*thrustAccelerationsPframe;


//    const Eigen::Vector3d thrustAccelerationsBframe = tudat::reference_frames::getPropulsionToBodyFrameTransformationMatrix(thrustAzimuthTest,thrustElevationTest)*thrustAccelerationsPframe;

    const Eigen::Vector3d thrustAccelerationsBframe = getPropulsionToBodyFrameTransformationMatrix(thrustAzimuthTest,thrustElevationTest)*thrustAccelerationsPframe;

/*
    std::cout<<"The thrustAccelerationsPframe is "<<thrustAccelerationsPframe<<std::endl;
    std::cout<<"The thrustAccelerationsBframe is "<<thrustAccelerationsBframe<<std::endl;
*/

    //// Testing the Auxiliary class ///

    /*
//    // Testing the atan2 function of c++

//    double atan2Test1 = std::atan2(0,0);        // Results in 0, but is in-fact undefined!!
//    double atan2Test2 = std::atan2(4,3);

//    std::cout<<"The atan2 of y = 0 and x = 0 is "<<atan2Test1<<std::endl;
//    std::cout<<"The atan2 of y = 4 and x = 3 is "<<atan2Test2<<std::endl;
*/


    Auxiliary Aux(adiabeticIndex, specificGasConstant,standardGravitationalParameter, rotationalVelocity, primeMeridianAngle,
              inertialFrameTime, bodyReferenceRadius, temperaturePolyCoefficients, temperatureAltitudeRanges,
              densityPolyCoefficients, Thrust, specificImpulse,
              referenceArea, dragCoefficientPolyCoefficients, dragCoefficientMachRanges);


    // Compute the auxiliary equations

    Eigen::VectorXd auxiliaryEquations =  Aux.getAuxiliaryEquations(aState,currentTime,thrustAccelerationsBframe);


    // Compute the auxiliary derivatives

    Eigen::VectorXd auxiliaryDerivatives = Aux.getAuxiliaryDerivatives(aState,currentTime,thrustAccelerationsBframe,auxiliaryEquations);

    // Compute the auxiliary functions

    Eigen::MatrixXd auxiliaryFunctions = Aux.getAuxiliaryFunctions(aState,currentTime,thrustAccelerationsBframe,auxiliaryEquations,auxiliaryDerivatives);

///*
    std::cout<<"The auxiliaryEquations are "<<auxiliaryEquations<<std::endl;
//    std::cout<<"The auxiliaryDerivatives are "<<auxiliaryDerivatives<<std::endl;
//    std::cout<<"The auxiliaryFunctions are "<<auxiliaryFunctions<<std::endl;
//*/

/*
    std::cout<<"The computed initial latitude is "<<auxiliaryEquations(12)<<std::endl;
    std::cout<<"The original initial latitude is "<<initialLatitude<<std::endl;
    std::cout<<"The computed initial latitude in deg is "<<rad2deg(auxiliaryEquations(12))<<std::endl;
    std::cout<<"The difference in rad is "<<initialLatitude-auxiliaryEquations(12)<<std::endl;
    std::cout<<"The difference in deg is "<<initialLatitudeDeg-rad2deg(auxiliaryEquations(12))<<std::endl;
 */
/*
    std::cout<<"Difference in radius = "<<initialRadius-auxiliaryEquations(20)<<std::endl;
    std::cout<<"initialRadius minus computed initial position inertial = "<<initialRadius-sqrt(aState(0)*aState(0)+aState(1)*aState(1)+aState(2)*aState(2))<<std::endl;
    std::cout<<"initialRadius minus computed initial position rotational = "<<initialRadius-sqrt(initialCartesianPositionRotationalFrame(0)*initialCartesianPositionRotationalFrame(0)+
                                                                                                 initialCartesianPositionRotationalFrame(1)*initialCartesianPositionRotationalFrame(1)+
                                                                                                 initialCartesianPositionRotationalFrame(2)*initialCartesianPositionRotationalFrame(2))<<std::endl;
*/

/*
    /// Debug ///

    std::cout<<"w24,1 = "<<auxiliaryFunctions(24,1)<<std::endl;
    std::cout<<"w24,2 = "<<auxiliaryFunctions(24,2)<<std::endl;
    std::cout<<"w24,3 = "<<auxiliaryFunctions(24,3)<<std::endl;
    std::cout<<"w24,4 = "<<auxiliaryFunctions(24,4)<<std::endl;
    std::cout<<"w24,5 = "<<auxiliaryFunctions(24,5)<<std::endl;
    std::cout<<"w24,6 = "<<auxiliaryFunctions(24,6)<<std::endl;
    std::cout<<"w24,7 = "<<auxiliaryFunctions(24,7)<<std::endl;
    std::cout<<"w24,8 = "<<auxiliaryFunctions(24,8)<<std::endl;
    std::cout<<"w24,9 = "<<auxiliaryFunctions(24,9)<<std::endl;
    std::cout<<"w24,10 = "<<auxiliaryFunctions(24,10)<<std::endl;
    std::cout<<"w24,11 = "<<auxiliaryFunctions(24,11)<<std::endl;
    std::cout<<"w24,12 = "<<auxiliaryFunctions(24,12)<<std::endl;
    std::cout<<"w24,13 = "<<auxiliaryFunctions(24,13)<<std::endl;
    std::cout<<"w24,14 = "<<auxiliaryFunctions(24,14)<<std::endl;
    std::cout<<"w24,15 = "<<auxiliaryFunctions(24,15)<<std::endl;
    std::cout<<"w24,16 = "<<auxiliaryFunctions(24,16)<<std::endl;
    std::cout<<"w24,17 = "<<auxiliaryFunctions(24,17)<<std::endl;

*/


    return 0;
}
