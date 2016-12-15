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
//    std::cout<<"The auxiliaryEquations are "<<auxiliaryEquations<<std::endl;
//    std::cout<<"The auxiliaryDerivatives are "<<auxiliaryDerivatives<<std::endl;
//    std::cout<<"The auxiliaryFunctions are "<<auxiliaryFunctions<<std::endl;
//*/

//    std::cout<<"This should be the size of the vector = "<<auxiliaryEquations.size()<<std::endl;

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




    //// Testing the basic recurrence relations ////

/*
    Eigen::VectorXd F = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd G = Eigen::VectorXd::Zero(4);

    int count = 0;

    for (int i = 0; i < F.size(); i++){

        F(i) = 1+i;
        G(i) = 2+i;

//        std::cout<<" F = "<<F<<std::endl;
//        std::cout<<" G = "<<G<<std::endl;

        count++;

//                std::cout<<"Count = "<<count<<std::endl;
    };


    /// Multiplication ///
    double Wmult = getMultiplicationRecurrenceRelation(F,G);


    std::cout<<"F test vector = "<<F<<std::endl;
    std::cout<<"G test vector = "<<G<<std::endl;
    std::cout<<"Wmult test value = "<<Wmult<<std::endl;




    /// Division ///
    Eigen::VectorXd WdivVector(1);              // Initialising the 0th order solution and vector

    WdivVector <<  F(0)/G(0);

//    std::cout<<"WdivVector = "<<WdivVector<<std::endl;



    for (int i=1; i<4; i++){

        double Wdiv = getDivisionRecurrenceRelation(F.segment(0,i+1),G.segment(0,i+1),WdivVector);      // Calling the recurrence relation using the first i+1 values for F and G and the vector of output data

//        std::cout<<"F.segment(0,"<<i+1<<") = "<<F.segment(0,i+1)<<std::endl;
//        std::cout<<"G.segment(0,"<<i+1<<") = "<<G.segment(0,i+1)<<std::endl;
//        std::cout<<"WdivVector = "<<WdivVector<<std::endl;

        Eigen::VectorXd intVect(i+1);               // Creating a bigger vector

        intVect << WdivVector.segment(0,i), Wdiv;       // Filling that bigger vector

        WdivVector = intVect;                       // Updating the output vector

//        std::cout<<"The output for WdivVector = "<<WdivVector<<std::endl;

    };

    std::cout<<"The output for WdivVector = "<<WdivVector<<std::endl;

    /// Power ///
    Eigen::VectorXd WpowVector(1);              // Initialising the 0th order solution and vector

    WpowVector <<  F(0)*F(0);                   // Perform the initial calculation (is already known through the w values)

    const double power = 2;                     // Set the power



    for (int i=1; i<4; i++){

        double Wpow = getPowerRecurrenceRelation(F.segment(0,i+1),WpowVector,power);      // Calling the recurrence relation using the first i+1 values for F and the vector of output data



        Eigen::VectorXd intVect(i+1);               // Creating a bigger vector

        intVect << WpowVector.segment(0,i), Wpow;       // Filling that bigger vector

        WpowVector = intVect;                       // Updating the output vector


    };

    std::cout<<"The output for WpowVector = "<<WpowVector<<std::endl;


    /// Exponential ///
    Eigen::VectorXd WexpVector(1);              // Initialising the 0th order solution and vector

    WexpVector <<  exp(F(0));                   // Perform the initial calculation (is already known through the w values)



    for (int i=1; i<4; i++){

        double Wexp = getExponentialRecurrenceRelation(F.segment(0,i+1),WexpVector);      // Calling the recurrence relation using the first i+1 values for F and the vector of output data



        Eigen::VectorXd intVect(i+1);               // Creating a bigger vector

        intVect << WexpVector.segment(0,i), Wexp;       // Filling that bigger vector

        WexpVector = intVect;                       // Updating the output vector


    };

    std::cout<<"The output for WexpVector = "<<WexpVector<<std::endl;







    /// Sine and Cosine ///

//             const double& pi = tudat::mathematical_constants::LONG_PI;


//            Eigen::VectorXd H(4);           // Creating a new testing vector with radians
//            H << pi/6, pi/3, pi/2, pi;      // Update: that didn't really help actually...

//            std::cout<<"H = "<<H<<std::endl;

    Eigen::VectorXd WcosVector(1);              // Initialising the 0th order solution and vector
    Eigen::VectorXd WsinVector(1);


    WcosVector <<  cos(F(0));                   // Perform the initial calculation (is already known through the w values)
    WsinVector <<  sin(F(0));


    for (int i=1; i<4; i++){

        double Wcos = getCosineRecurrenceRelation(F.segment(0,i+1),WsinVector);      // Calling the recurrence relation using the first i+1 values for F and the vector of output data
        double Wsin = getSineRecurrenceRelation(F.segment(0,i+1),WcosVector);


        Eigen::VectorXd intVectCos(i+1);               // Creating a bigger vector
        Eigen::VectorXd intVectSin(i+1);

        intVectCos << WcosVector.segment(0,i), Wcos;       // Filling that bigger vector
        intVectSin << WsinVector.segment(0,i), Wsin;

        WcosVector = intVectCos;                       // Updating the output vector
        WsinVector = intVectSin;


    };

    std::cout<<"The output for WcosVector = "<<WcosVector<<std::endl;
    std::cout<<"The output for WsinVector = "<<WsinVector<<std::endl;

    //*/

    /// Testing all the recurrence relations ///

    const int maxOrder = 20;
///*
    Eigen::MatrixXd TaylorCoefficients = getTaylorCoefficients(adiabeticIndex, specificGasConstant, standardGravitationalParameter, rotationalVelocity, primeMeridianAngle,
                          inertialFrameTime, bodyReferenceRadius,temperaturePolyCoefficients, temperatureAltitudeRanges,
                          densityPolyCoefficients, Thrust, specificImpulse,
                          referenceArea, dragCoefficientPolyCoefficients, dragCoefficientMachRanges,
            thrustAccelerationsBframe,
            auxiliaryEquations,
            auxiliaryDerivatives,
            auxiliaryFunctions,
            currentTime,
            maxOrder);

//        std::cout<<"The Taylor Coefficients are: "<<TaylorCoefficients<<std::endl;

//*/

    /// Testing the write capability to store the Taylor Coefficients ///




        Eigen::MatrixXd TaylorCoefficientsOutputMatrix = Eigen::MatrixXd::Zero(7,maxOrder+1);       // Create an output matrix for the file without the first empty row

        TaylorCoefficientsOutputMatrix.row(0) = TaylorCoefficients.row(1);                  // The first line entries are the maxOrder+1 Taylor Series Coefficients for     the position in the x-direction
        TaylorCoefficientsOutputMatrix.row(1) = TaylorCoefficients.row(2);                  // The second line entries are the maxOrder+1 Taylor Series Coefficients for    the position in the y-direction
        TaylorCoefficientsOutputMatrix.row(2) = TaylorCoefficients.row(3);                  // The third line entries are the maxOrder+1 Taylor Series Coefficients for     the position in the z-direction
        TaylorCoefficientsOutputMatrix.row(3) = TaylorCoefficients.row(4);                  // The fourth line entries are the maxOrder+1 Taylor Series Coefficients for    the velocity in the x-direction
        TaylorCoefficientsOutputMatrix.row(4) = TaylorCoefficients.row(5);                  // The fifth line entries are the maxOrder+1 Taylor Series Coefficients for     the velocity in the y-direction
        TaylorCoefficientsOutputMatrix.row(5) = TaylorCoefficients.row(6);                  // The sixth line entries are the maxOrder+1 Taylor Series Coefficients for     the velocity in the z-direction
        TaylorCoefficientsOutputMatrix.row(6) = TaylorCoefficients.row(7);                  // The seventh line entries are the maxOrder+1 Taylor Series Coefficients for   the mass



        // Set directory where output files will be stored. By default, this is your project
        // root-directory.
        const std::string outputDirectory = "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/tudatApplications/thesisProject/testOutputFolder/";

//        std::cout<<"The output directory = "<<outputDirectory<<std::endl;



        // Set output format for matrix output.
        Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

        // Set absolute path to file containing the Taylor Series Coefficients.
        const std::string taylorSeriesCoefficientsAbsolutePath = outputDirectory + "test4TaylorSeriesCoefficients.csv";


        // Check if the file already exists.


        std::ifstream ifile(taylorSeriesCoefficientsAbsolutePath.c_str()); // Check it as an input file

        bool fexists = false;

        if (ifile){


           fexists = true;

           ifile.close();

        }


//        std::cout<<"Does the file exist? (0 = no, 1 = yes) "<<fexists<<std::endl;



        // If so: append, if not: create new file and put data in

        if (fexists == true){

            // Export the Taylor Series Coefficients matrix.
            std::ofstream exportFile1;                          // Define the file as an output file


            exportFile1.open(taylorSeriesCoefficientsAbsolutePath.c_str(),std::ios_base::app);      // Open the file in append mode

            exportFile1 << "\n";                                            // Make sure the new matrix

            exportFile1 << TaylorCoefficientsOutputMatrix.format( csvFormat );


            exportFile1.close( );
}
            else{

            // Export the Taylor Series Coefficients matrix.
            std::ofstream exportFile1( taylorSeriesCoefficientsAbsolutePath.c_str( ) );
            exportFile1 << TaylorCoefficientsOutputMatrix.format( csvFormat );
            exportFile1.close( );
        };


        /// Testing the StepSize class ///

/*
        // Test to see how the order of magnitude works
        double de = 1.2345E-5;
        int orderOfMagnitude = floor(log10(de));

        std::cout<<"The order of magnitude of de = "<<orderOfMagnitude<<std::endl;

//*/
        tudat::basic_mathematics::Vector7d penultimateCoefficients;
        tudat::basic_mathematics::Vector7d lastCoefficients;

        for (int i = 0; i<7; i++){

            penultimateCoefficients(i)= TaylorCoefficients((i+1),maxOrder-1);
            lastCoefficients(i)= TaylorCoefficients((i+1),(maxOrder));
        }

//        std::cout<<"The penultimateCoefficients are "<<penultimateCoefficients<<std::endl;
//        std::cout<<"The lastCoefficients are "<<lastCoefficients<<std::endl;


        StepSize stepSize; // Initializing the stepSize class. THIS SHOULD BE DONE BEFORE THE START OF THE INTEGRATION!!!!!

        // Checking the default values

        double currentStepSize = stepSize.getCurrentStepSize();
        double localErrorTolerance = stepSize.getLocalErrorTolerance();
        double stepMultiplicationFactor = stepSize.getStepMultiplicationFactor();

        std::cout<<"The current step-size = "<<currentStepSize<<std::endl;
        std::cout<<"The local error tolerance = "<<localErrorTolerance<<std::endl;
        std::cout<<"The step multiplication factor = "<<stepMultiplicationFactor<<std::endl;
//*/

 /*       // Adjusting the values

        stepSize.setCurrentStepSize(2.3);
        stepSize.setLocalErrorTolerance(1E-7);
        stepSize.setStepMultiplicationFactor(0.8);

        currentStepSize = stepSize.getCurrentStepSize();
        localErrorTolerance = stepSize.getLocalErrorTolerance();
        stepMultiplicationFactor = stepSize.getStepMultiplicationFactor();

        std::cout<<"The updated step-size = "<<currentStepSize<<std::endl;
        std::cout<<"The updated local error tolerance = "<<localErrorTolerance<<std::endl;
        std::cout<<"The updated step multiplication factor = "<<stepMultiplicationFactor<<std::endl;

        //*/

/*        // Updating the step-size using the previous step-size
        stepSize.updateStepSizeUsingPreviousStepSize(penultimateCoefficients, lastCoefficients, maxOrder);

        currentStepSize = stepSize.getCurrentStepSize();
        localErrorTolerance = stepSize.getLocalErrorTolerance();
        stepMultiplicationFactor = stepSize.getStepMultiplicationFactor();
        double maximumTruncationErrorEstimate = stepSize.getMaximumTruncationErrorEstimate();
//        double orderMaxTruncErrorEstimate = stepSize.getOrderMaxTruncErrorEstimate();
//        tudat::basic_mathematics::Vector7d truncationErrorEstimates = stepSize.getTruncationErrorEstimates();

        std::cout<<"The current step-size from previous step-size method = "<<currentStepSize<<std::endl;
//        std::cout<<"The local error tolerance = "<<localErrorTolerance<<std::endl;
//        std::cout<<"The step multiplication factor = "<<stepMultiplicationFactor<<std::endl;

        std::cout<<"The maximumTruncationErrorEstimate = "<<maximumTruncationErrorEstimate<<std::endl;
//        std::cout<<"The order of the maximumTruncationErrorEstimate = "<<orderMaxTruncErrorEstimate<<std::endl;
//        std::cout<<"The truncation error estimates = "<<truncationErrorEstimates<<std::endl;
//*/


        /// Testing the actual Taylor Series expansion for every state variable ///

        tudat::basic_mathematics::Vector7d updatedState;        // Create a vector for the updatedState

        for (int n = 0; n<updatedState.size();n++){                 // All variables

        for (int k = 0; k<maxOrder+1;k++){                      // Taylor series summation

            updatedState(n) += TaylorCoefficients((n+1),k)*pow(currentStepSize,k);      // Perform one step of the taylor series expansion and then add it to the previous step

        } // Taylor series summation

}   // All variables



        double updatedTime = currentTime+currentStepSize;           // Create the updated time variable

        std::cout<<"updatedTime = "<<updatedTime<<std::endl;

        std::cout<<"updatedState = "<<updatedState<<std::endl;




///*        // Updating the step-size using the iteration method (This method was perferred by Scott and Martini [2008]


       stepSize.updateStepSizeUsingIteration(penultimateCoefficients, lastCoefficients, maxOrder);

       currentStepSize = stepSize.getCurrentStepSize();
       localErrorTolerance = stepSize.getLocalErrorTolerance();
       stepMultiplicationFactor = stepSize.getStepMultiplicationFactor();



       std::cout<<"The current step-size from iteration method = "<<currentStepSize<<std::endl;
       std::cout<<"The local error tolerance = "<<localErrorTolerance<<std::endl;
       std::cout<<"The step multiplication factor = "<<stepMultiplicationFactor<<std::endl;
//*/


    return 0;
}


