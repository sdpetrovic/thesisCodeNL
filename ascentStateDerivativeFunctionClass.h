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
 *      160518    S.D. Petrovic     File created
 *      160705    S.D. Petrovic     Added rounding error corrections
 *
 *    References
 *
 *    Notes
 *
 */


#ifndef ASCENTSTATEDERIVATIVEFUNCTIONCLASS_H
#define ASCENTSTATEDERIVATIVEFUNCTIONCLASS_H


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
#include <thesisProject/celestialBody.h>            // Final version

/// Testing the vehicle class ///
#include <thesisProject/MarsAscentVehicle.h>    // Final version

/// Testing the current state and time and its updater ///
#include <thesisProject/stateAndTime.h>             // Final version

/// Testing the auxiliary equations ///
//#include <thesisProject/Auxiliary.h>                // Original test file

/// Testing the other required functions ///
#include <thesisProject/projectLibraries/otherRequiredFunctions.h>              // Original test file

class ascentStateDerivativeFunctionClass
{
public:
    ascentStateDerivativeFunctionClass(const celestialBody& planet_, const MarsAscentVehicle& MAV_){

        // Create variables to be used in this function

            Mars = planet_;                               // The celestial body class
            MAV = MAV_;                                 // The MAV class
//            StateAndTime stateAndTime = currentStateAndTime_;      // The current state and time class

            rotationalVelocityMars = Mars.rotationalVelocity();
            primeMeridianAngle = Mars.primeMeridianAngle();
            inertialFrameTime = Mars.inertialFrameTime();

            thrustAzimuth = MAV.thrustAzimuth();
            thrustElevation = MAV.thrustElevation();


            FlightPathAngle = 1000;
            HeadingAngle = 1000;


    }         // NOTE TO SELF: WHEN INITIALIZING A NEW CLASS WITHOUT ANY PRESETS, MAKE SURE TO INCLUDE {} BEHIND (). IN THIS CASE THAT WOULD LOOK LIKE THIS: ascentStateDerivativeFunctionClass(){}

    /// Debug ///

    const double getTestValue(){

        const double testValue = 1;

        return testValue;
    }

    /// Debug ///

    /// Set functions ///

            void setFlightPathAngleAndHeadingAngle(const double FlightPathAngle_ = 1000, const double HeadingAngle_ = 1000){
                if (FlightPathAngle_ >= -tudat::mathematical_constants::LONG_PI/2.0 && FlightPathAngle_ <= tudat::mathematical_constants::LONG_PI/2.0){
                    if (HeadingAngle_ >= 0.0 && HeadingAngle_ <= 2.0*tudat::mathematical_constants::LONG_PI){
                        FlightPathAngle = FlightPathAngle_; // Flight path angle in rad
                        HeadingAngle = HeadingAngle_;  // Heading angle in rad

                        std::cout<<"Flight-path angle and heading angle have been set"<<std::endl;
                    }
                    else{
                        std::cout<<"Heading angle has to be specified between 0.0 and 2*pi"<<std::endl;
                    }

                }
                else if (HeadingAngle_ >= 0.0 && HeadingAngle_ <= 2.0*tudat::mathematical_constants::LONG_PI) {
                    std::cout<<"Flight-path angle has to be specified between -pi/2 and pi/2"<<std::endl;
                }
                else {
                    std::cout<<"Heading angle has to be specified between 0.0 and 2*pi and Flight-path angle has to be specified between -pi/2 and pi/2"<<std::endl;
                }

            }



    tudat::basic_mathematics::Vector7d ascentStateDerivativeFunction(const double currentTime_, const tudat::basic_mathematics::Vector7d& currentState){

        /// Debug ///
        //std::cout<<"Here 0"<<std::endl;
        /// Debug ///

    /// Compute current spherical state ///

//        tudat::basic_mathematics::Vector7d currentState = stateAndTime.getCurrentState();

        const double xPosition = currentState(0);            // x position coordinate definition
        const double yPosition = currentState(1);            // y position coordinate definition
        const double zPosition = currentState(2);            // z position coordinate definition
        const double xVelocity = currentState(3);            // x velocity coordinate definition
        const double yVelocity = currentState(4);            // y velocity coordinate definition
        const double zVelocity = currentState(5);            // z velocity coordinate definition
        const double massMAV = currentState(6);              // MAV mass definition
        const double currentTime = currentTime_;   // current time definition

    // Computations

        const double Radius = sqrt(xPosition*xPosition+yPosition*yPosition+zPosition*zPosition);         // r [km]

        const double inertialVelocity = sqrt(xVelocity*xVelocity+yVelocity*yVelocity+zVelocity*zVelocity);       // V_I [km/s]

//        const double inertialLongitude = atan2(yPosition,xPosition);         // lambda [rad]

        const double Latitude = asin(zPosition/Radius);              // delta [rad]

        /// Debug ///
        //std::cout<<"Here 1"<<std::endl;
        /// Debug ///

//        const double rotationalLongitude = inertialLongitude-rotationalVelocityMars*(inertialFrameTime+currentTime)+primeMeridianAngle;  // tau [rad]

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


        double rotationalVelocity_;   // V_G [km/s] (placeholder)
        // Deal with rounding errors
        if (inertialVelocity*inertialVelocity+rotationalVelocityMars*rotationalVelocityMars*(xPosition*xPosition+yPosition*yPosition)+2.0*rotationalVelocityMars*(xVelocity*yPosition-yVelocity*xPosition) <= 0.0){
            rotationalVelocity_ = 0.0;
        }
        else {
         rotationalVelocity_ = sqrt(inertialVelocity*inertialVelocity+rotationalVelocityMars*rotationalVelocityMars*(xPosition*xPosition+yPosition*yPosition)+2.0*rotationalVelocityMars*(xVelocity*yPosition-yVelocity*xPosition)); // V_G [km/s]
        }

        const double rotationalVelocity = rotationalVelocity_;      // V_G [km/s] (actual parameter)

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
        /// Debug ///
        //std::cout<<"Here 4"<<std::endl;
        /// Debug ///

        double rotationalAzimuth_;
        double rotationalFlightPathAngle_;          // gamma_G [rad] (placeholder)

        if (FlightPathAngle != 1000 && HeadingAngle != 1000 && currentTime == 0.0){
            rotationalAzimuth_ = HeadingAngle;
            rotationalFlightPathAngle_ = FlightPathAngle;

        }
        else {

            rotationalAzimuth_ = atan2(verticalYvelocity,verticalXvelocity);    // chi_G [rad]

            // Avoid singularities
            if (rotationalVelocity == 0.0){
                rotationalFlightPathAngle_ = tudat::mathematical_constants::LONG_PI/2.0;
            }
            else if (verticalZvelocity/rotationalVelocity >= 1.0 || verticalZvelocity/rotationalVelocity-1 >= -1E-15){  // Compensate for rounding errors
    //            std::cout<<"sin(FPA) has been rounded down to 1 with difference: "<<abs(verticalZvelocity/rotationalVelocity)-1<<std::endl;
                rotationalFlightPathAngle_ = -asin(1.0);
            }
            else if (verticalZvelocity/rotationalVelocity <= -1.0 || verticalZvelocity/rotationalVelocity+1 <= 1E-15){ // Compensate for rounding errors
                rotationalFlightPathAngle_ = -asin(-1.0);
            }
            else {
            rotationalFlightPathAngle_ = -asin(verticalZvelocity/rotationalVelocity);   // gamma_G [rad]

            }
        }

//        const double rotationalAzimuth = atan2(verticalYvelocity,verticalXvelocity);    // chi_G [rad]

        const double rotationalAzimuth = rotationalAzimuth_;        // chi_G [rad]


        const double rotationalFlightPathAngle = rotationalFlightPathAngle_; // gamma_G [rad] (actual parameter)



        /// Debug ///
//        std::cout<<"Here 5"<<std::endl;
//        std::cout<<"Latitude = "<<Latitude<<std::endl;
//        std::cout<<"rotationalLongitude = "<<rotationalLongitude<<std::endl;
//        std::cout<<"FlightPathAngle = "<<rotationalFlightPathAngle<<std::endl;
//        std::cout<<"FlightPathAngle (deg) = "<<rad2deg(rotationalFlightPathAngle)<<std::endl;
//        std::cout<<"HeadingAngle = "<<rotationalAzimuth<<std::endl;
//        std::cout<<"HeadingAngle (deg) = "<<rad2deg(rotationalAzimuth)<<std::endl;

//        std::cout<<"RotationalVelocity = "<<rotationalVelocity<<std::endl;
//        std::cout<<"verticalXvelocity = "<<verticalXvelocity<<std::endl;
//        std::cout<<"verticalYvelocity = "<<verticalYvelocity<<std::endl;
//        std::cout<<"verticalZvelocity = "<<verticalZvelocity<<std::endl;
//        std::cout<<"rotationalVelocity = "<<rotationalVelocity<<std::endl;
//        std::cout<<"Radius = "<<Radius<<std::endl;



        /// Debug ///

        /* // Old functions
//        double inertialLongitudeChange_;       // lambda_dot [rad/s] (placeholder)

//        // Avoiding singularities
//        if ((xPosition*xPosition+yPosition*yPosition) == 0){

//            inertialLongitudeChange_ = 0.0;
//        }
//        else {
//            inertialLongitudeChange_ = (xPosition*yVelocity-yPosition*xVelocity)/(xPosition*xPosition+yPosition*yPosition);
//        };

//        const double inertialLongitudeChange = inertialLongitudeChange_;        // lambda_dot [rad/s] (actual parameter)


//        double rotationalLongitudeChange_ = inertialLongitudeChange-rotationalVelocityMars;     // tau_dot [rad/s] (placeholder)

//        if (abs(rotationalLongitudeChange_)<=1e-11){  // Setting the accuracy to 1e-15 to avoid problems in the beginning with rounding errors...

//            rotationalLongitudeChange_ = 0.0;

//        };

//        const double rotationalLongitudeChange = rotationalLongitudeChange_;    // tau_dot [rad/s] (actual parameter)

        /// Debug ///

        std::cout<<"rotationalVelocityMars = "<<rotationalVelocityMars<<std::endl;
        std::cout<<"inertialLongitudeChange = "<<inertialLongitudeChange<<std::endl;
        std::cout<<"rotationalVelocityMars-7.088e-05 = "<<rotationalVelocityMars-7.088e-05<<std::endl;
        std::cout<<"inertialLongitudeChange-7.088e-05 = "<<inertialLongitudeChange-7.088e-05<<std::endl;
        std::cout<<"(xPosition*yVelocity-yPosition*xVelocity) = "<<(xPosition*yVelocity-yPosition*xVelocity)<<std::endl;
        std::cout<<"(xPosition*xPosition+yPosition*yPosition) = "<<(xPosition*xPosition+yPosition*yPosition)<<std::endl;



//        const double RadiusChange = (xPosition*xVelocity+yPosition*yVelocity+zPosition*zVelocity)/(Radius);      // radial velocity [km/s]

//        double LatitudeChange_; // delta_dot [rad/s] (placeholder)

//        if ((Radius*Radius*sqrt(1-(zPosition/Radius)*(zPosition/Radius))) == 0){
//            LatitudeChange_ = 0.0;
//        }
//        else{
//            LatitudeChange_ = (Radius*zVelocity-zPosition*RadiusChange)/(Radius*Radius*sqrt(1.0-(zPosition/Radius)*(zPosition/Radius)));
//        };

//        const double LatitudeChange = LatitudeChange_;  // delta_dot [rad/s] (actual parameter)

//        double localMarsRotationalVelocity = rotationalVelocityMars*Radius*cx12;  // V_M [km/s]
//        // Avoid cosine rounding errors
//        if (abs(cx12)<6.2e-17){
//          localMarsRotationalVelocity = 0.0;
//        }

//        double inertialFlightPathAngle = asin(RadiusChange/inertialVelocity);          // gamma_I [rad]
//        // Avoid singularities
//        if (inertialVelocity == 0){
//            inertialFlightPathAngle = 0.0;
//        }
//        // And dealing with round-off errors
//        else if ((RadiusChange/inertialVelocity) < -1.0 && ((RadiusChange/inertialVelocity)+1.0) > -1e-15){
//            inertialFlightPathAngle = asin(-1.0);
////            std::cout<<"Rounding to -1"<<std::endl;
//        }
//        else if ((RadiusChange/inertialVelocity) > 1.0 && ((RadiusChange/inertialVelocity)-1.0) < 1e-15){
//            inertialFlightPathAngle = asin(1.0);
////            std::cout<<"Rounding to 1"<<std::endl;
//        }

//        const double inertialAzimuth = atan2((inertialLongitudeChange*cx12),LatitudeChange);    // chi_I [rad]

//        const double rotationalVelocity = sqrt(localMarsRotationalVelocity*localMarsRotationalVelocity+inertialVelocity*inertialVelocity-2.0*localMarsRotationalVelocity*inertialVelocity*cos(inertialFlightPathAngle)*sin(inertialAzimuth));  // V_R [km/s]

//        double rotationalFlightPathAngle_; // gamma_R [rad]  (placeholder)

//        if (rotationalVelocity == 0){       // Setting the initial flight path angle in the rotational frame to 90 deg (or pi/s)

//            rotationalFlightPathAngle_ = tudat::mathematical_constants::LONG_PI/2.0;
////            rotationalFlightPathAngle_ = pi/2.0;
//        }
//        else if ((RadiusChange/rotationalVelocity)<-1.0 && ((RadiusChange/rotationalVelocity)+1.0)>-1e-15){
//            rotationalFlightPathAngle_=asin(-1.0);
//        }
//        else if ((RadiusChange/rotationalVelocity)>1.0 && ((RadiusChange/rotationalVelocity)-1.0)<1e-15){
//            rotationalFlightPathAngle_ = asin(1.0);
//        }
//        else {
//            rotationalFlightPathAngle_ = asin(RadiusChange/rotationalVelocity);
//        };

//        const double rotationalFlightPathAngle = rotationalFlightPathAngle_;    // gamma_R [rad] (actual parameter)



//        const double rotationalAzimuth = atan2((rotationalLongitudeChange*cx12),LatitudeChange);    // chi_R [rad]


       // Check output
        std::cout<<"///////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
        std::cout<<"Radius = "<<Radius<<std::endl;
        std::cout<<"inertialVelocity = "<<inertialVelocity<<std::endl;
//        std::cout<<"inertialLongitude = "<<inertialLongitude<<std::endl;
        std::cout<<"Latitude = "<<Latitude<<std::endl;
        std::cout<<"rotationalLongitude = "<<rotationalLongitude<<std::endl;
//        std::cout<<"inertialLongitudeChange = "<<inertialLongitudeChange<<std::endl;
//        std::cout<<"rotationalLongitudeChange = "<<rotationalLongitudeChange<<std::endl;
//        std::cout<<"RadiusChange = "<<RadiusChange<<std::endl;
//        std::cout<<"RadiusChange/rotationalVelocity = "<<RadiusChange/rotationalVelocity<<std::endl;
//        std::cout<<"RadiusChange/rotationalVelocity+1 = "<<RadiusChange/rotationalVelocity+1<<std::endl;
//        std::cout<<"asin(RadiusChange/rotationalVelocity) = "<<asin(RadiusChange/rotationalVelocity)<<std::endl;
//        std::cout<<"LatitudeChange = "<<LatitudeChange<<std::endl;
//        std::cout<<"localMarsRotationalVelocity = "<<localMarsRotationalVelocity<<std::endl;
//        std::cout<<"inertialFlightPathAngle = "<<inertialFlightPathAngle<<std::endl;
//        std::cout<<"inertialAzimuth = "<<inertialAzimuth<<std::endl;
        std::cout<<"rotationalVelocity = "<<rotationalVelocity<<std::endl;
        std::cout<<"rotationalFlightPathAngle = "<<rotationalFlightPathAngle<<std::endl;
        std::cout<<"rotationalAzimuth = "<<rotationalAzimuth<<std::endl;
        std::cout<<"///////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
//*/



    /// Testing the local air temperature function ///

        const double currentAltitude = Radius-Mars.bodyReferenceRadius();

        const double currentTemperature = air_temperature::airTemperature(Mars.temperaturePolyCoefficients(), Mars.temperatureAltitudeRanges(),currentAltitude);

 /*      // Check output
        std::cout<<"currentAltitude = "<<currentAltitude<<std::endl;
        std::cout<<"currentTemperature = "<<currentTemperature<<std::endl;
        std::cout<<"Radius = "<<Radius<<std::endl;
        std::cout<<"R_MOLA = "<<Mars.bodyReferenceRadius()<<std::endl;
        std::cout<<"Radius-3395.4 = "<<Radius-3395.4<<std::endl;
        std::cout<<"R_MOLA-3396 = "<<Mars.bodyReferenceRadius()-3396<<std::endl;
    //*/

    /// Testing the local air density function ///

        const double currentDensity= air_density::airDensity(Mars.densityPolyCoefficients(),  currentAltitude);

//        std::cout<<"The current air density = "<<currentDensity<<std::endl;

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

//        std::cout<<"currentDrag = "<<currentDrag<<std::endl;



        /// Thrust acceleration in B-frame ///   thrustAccelerationsBframe

        const Eigen::Vector3d thrustAccelerationsPframe = Eigen::Vector3d((MAV.Thrust()/massMAV),0,0);            // THIS HAS TO BE CHANGED IN THE FUTURE TO INCLUDE A WIDE RANGE OF THRUST AZIMUTH AND ELEVATION ANGLES!!!

//        const double thrustAzimuthTestDeg = 0.0;             // thrust azimuth gimbal angle [Deg] 10 for testing
//        const double thrustElevationTestDeg = 0.0;            // thrust elevation gimbal angle [Deg] 5 for testing

//        const double thrustAzimuthTest = deg2rad(thrustAzimuthTestDeg);     // thrust azimuth gimbal angle [rad]
//        const double thrustElevationTest = deg2rad(thrustElevationTestDeg); // thrust elevation gimbal angle [rad]



        // Determine the proper azimuth value for the current altitude section
        int sectionThrustAz = 0;    // Set the current azimuth value to the default first section
        for (int i = 0; i < thrustAzimuth.rows();i++){
            if (thrustAzimuth(i,0) <= currentAltitude && currentAltitude < thrustAzimuth(i,1)){ // Test for all the sections (independent of how many sections there are)
                sectionThrustAz = i;
//                std::cout<<"sectionThrustAz = "<<sectionThrustAz+1<<std::endl;
            }
        }

        const double thrustAzimuthTest = thrustAzimuth(sectionThrustAz,2); // Set the thrust azimuth to the current azimuth corresponding to the current altitude section

        // Determine the proper elevation value for the current altitude section
        int sectionThrustEl = 0;    // Set the current elevation value to the default first section
        for (int i = 0; i < thrustElevation.rows();i++){
            if (thrustElevation(i,0) <= currentAltitude && currentAltitude < thrustElevation(i,1)){ // Test for all the sections (independent of how many sections there are)
                sectionThrustEl = i;
            }
        }

        const double thrustElevationTest = thrustElevation(sectionThrustEl,2); // Set the thrust elevation to the current elevation corresponding to the current altitude section


        /// Debug ///
//        std::cout<<"/// Debug ///"<<std::endl;
//        std::cout<<"thrustAzimuth.cols() = "<<thrustAzimuth.cols()<<std::endl;
//        std::cout<<"thrustAzimuth.rows() = "<<thrustAzimuth.rows()<<std::endl;
//        std::cout<<"currentAltitude = "<<currentAltitude<<std::endl;
//        std::cout<<"thrustAzimuthTest = "<<thrustAzimuthTest<<std::endl;
//        std::cout<<"thrustElevationTest = "<<thrustElevationTest<<std::endl;

//        std::cout<<"/// Debug ///"<<std::endl;
//        std::cout<<" "<<std::endl;
        /// Debug ///



//        const Eigen::Vector3d thrustAccelerationsBframe = getPropulsionToBodyFrameTransformationMatrix(thrustAzimuthTest,thrustElevationTest)*thrustAccelerationsPframe;

        // Test
        Eigen::Vector3d thrustAccelerationsBframe = getPropulsionToBodyFrameTransformationMatrix(thrustAzimuthTest,thrustElevationTest)*thrustAccelerationsPframe;

//        thrustAccelerationsBframe(0) = 0.0;
//        thrustAccelerationsBframe(1) = 0.0;
//        thrustAccelerationsBframe(2) = 0.0;
       // Test


//        std::cout<<"The thrust accelerations in the B-frame are "<<thrustAccelerationsBframe<<std::endl;

        /// Drag acceleration in B-frame ///

        const Eigen::Vector3d dragAccelerationsBframe = Eigen::Vector3d((-currentDrag/massMAV),0,0);

//        std::cout<<"The drag accelerations in the B-frame are "<<dragAccelerationsBframe<<std::endl;

    /// Transfer to the inertial frame ///
    ///
        /// Accelerations in the V-frame ///

        // Thrust accelerations from B-frame to V-frame
        const Eigen::Vector3d thrustAccelerationsVframe = tudat::reference_frames::getTrajectoryToLocalVerticalFrameTransformationMatrix(rotationalFlightPathAngle,rotationalAzimuth)*thrustAccelerationsBframe;

        // Drag accelerations from B-frame to V-frame
        const Eigen::Vector3d dragAccelerationsVframe = tudat::reference_frames::getTrajectoryToLocalVerticalFrameTransformationMatrix(rotationalFlightPathAngle,rotationalAzimuth)*dragAccelerationsBframe;

//        std::cout<<"The thrust accelerations in the V-frame are "<<thrustAccelerationsVframe<<std::endl;
    //    std::cout<<"The drag accelerations in the V-frame are "<<dragAccelerationsVframe<<std::endl;


        /// Accelerations in the R-frame ///

        // Thrust accelerations from V-frame to R-frame
        const Eigen::Vector3d thrustAccelerationsRframe = tudat::reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(rotationalLongitude,Latitude)*thrustAccelerationsVframe;

        // Drag accelerations from V-frame to R-frame
        const Eigen::Vector3d dragAccelerationsRframe = tudat::reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(rotationalLongitude,Latitude)*dragAccelerationsVframe;

//        std::cout<<"The thrust accelerations in the R-frame are "<<thrustAccelerationsRframe<<std::endl;
    //    std::cout<<"The drag accelerations in the R-frame are "<<dragAccelerationsRframe<<std::endl;

        /// Accelerations in the I-frame ///

        const double angleItoR = rotationalVelocityMars*(inertialFrameTime+currentTime)-primeMeridianAngle;

        // Thrust accelerations from R-frame to I-frame
        const Eigen::Vector3d thrustAccelerationsIframe = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(angleItoR)*thrustAccelerationsRframe;

        // Drag accelerations from R-frame to I-frame
        const Eigen::Vector3d dragAccelerationsIframe = tudat::reference_frames::getRotatingPlanetocentricToInertialFrameTransformationMatrix(angleItoR)*dragAccelerationsRframe;

//        std::cout<<"The thrust accelerations in the I-frame are "<<thrustAccelerationsIframe<<std::endl;
//        std::cout<<thrustAccelerationsIframe(0)<<", "<<thrustAccelerationsIframe(1)<<", "<<thrustAccelerationsIframe(2)<<", "<<std::endl;
    //    std::cout<<"The drag accelerations in the I-frame are "<<dragAccelerationsIframe<<std::endl;

    //    std::cout<<"The thrust+drag accelerations in the I-frame are "<<thrustAccelerationsIframe+dragAccelerationsIframe<<std::endl;


    /// Compute gravitational acceleration ///

        Eigen::Vector3d gravAccelerationsIframe_;        // Define the placeholder gravitational acceleration vector

        for (int i = 0; i<3;i++){

            gravAccelerationsIframe_(i)=-Mars.standardGravitationalParameter()*(currentState(i)/(pow(Radius,3.0)));

    //        std::cout<<"gravAccelerationsIframe_("<<i<<") = "<<gravAccelerationsIframe_(i)<<std::endl;
        };

        const Eigen::Vector3d gravAccelerationsIframe = gravAccelerationsIframe_;         // Actual gravitational acceleration vector


//        std::cout<<"The gravitational accelerations in the I-frame are "<<gravAccelerationsIframe<<std::endl;

    /// Compute total acceleration ///

        const Eigen::Vector3d totalAccelerationsIframe = gravAccelerationsIframe+dragAccelerationsIframe+thrustAccelerationsIframe;



    /// Compute the mass flow rate ///

        const double massFlowRate = -MAV.Thrust()/(tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION*MAV.specificImpulse());



//        std::cout<<"The velocities in the I-frame are "<<currentState(3)<<"\n"<<
//                currentState(4)<<"\n"<<
//                currentState(5)<<std::endl;
//        std::cout<<"The total accelerations in the I-frame are "<<totalAccelerationsIframe<<std::endl;
//        std::cout<<"The mass flow rate = "<<massFlowRate<<std::endl;


    /// Define the output vector and fill it ///


        tudat::basic_mathematics::Vector7d stateDerivativeVector;

        stateDerivativeVector(0) = xVelocity;
        stateDerivativeVector(1) = yVelocity;
        stateDerivativeVector(2) = zVelocity;
        stateDerivativeVector(3) = totalAccelerationsIframe(0);
        stateDerivativeVector(4) = totalAccelerationsIframe(1);
        stateDerivativeVector(5) = totalAccelerationsIframe(2);
        stateDerivativeVector(6) = massFlowRate;

        /// Debug ///
//        std::cout<<"stateDerivativeVector = "<<stateDerivativeVector<<std::endl;
        /// Debug ///

        return stateDerivativeVector;
        } // end of the state derivative function


private:

    // Create variables to be used in this function

        celestialBody Mars;                               // The celestial body class       // This initilizes a new Mars, but later replaces it by the provided planet
        MarsAscentVehicle MAV;                                 // The MAV class             // This initializes a new MAV, but later replaces it by the provided MAV

        double rotationalVelocityMars;                  // The rotational velocity of Mars [rad/s]
        double primeMeridianAngle;                      // The angle of the prime Meridian of Mars at time of inertial frame set [rad]
        double inertialFrameTime;                       // The time at inertial frame set [s]

        Eigen::MatrixXd thrustAzimuth;                  // The thrust azimuth angle values and the altitude boundary values [rad] and [km]
        Eigen::MatrixXd thrustElevation;                // The thrust elevation angle values and the altitude boundary values [rad] and [km]


        // Set functions

      double FlightPathAngle;         // Flight path angle in rad
     double HeadingAngle;            // Heading angle in rad


}; // end of class

#endif // ASCENTSTATEDERIVATIVEFUNCTIONCLASS_H
