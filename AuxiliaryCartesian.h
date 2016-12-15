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
 *      160413    S.D. Petrovic     File created
 *      160518    S.D. Petrovic     Fixed the mistake I made with the transformation matrix T_IB which resulted in a mistake in u6
 *      160520    S.D. Petrovic     Fixed mistake in x25 where it said x25 = x25/(2*x20) which should be x25 = x26/(2*x20). Also, u41 had a - instead of a +.
 *                                  Also, at u45 and u21 the tolerance had to include an abs function!
 *      160526    S.D. Petrovic     Added W4,0 to be able to properly evaluate the recurrence relation of W4,2
 *      160527    S.D. Petrovic     Corrected mistake in x42
 *      160531    S.D. Petrovic     Added more ways to deal with rounding errors
 *      160602    S.D. Petrovic     Added the thrust auxiliary functions
 *      160603    S.D. Petrovic     Updated u10 to be OmegaM instead of 0!
 *      160618    S.D. Petrovic     Found mistake in u24 and updated it and the corresponding auxiliary functions
 *      160622    S.D. Petrovic     Found a huge mistake in u15 (extra acc were not taken into account, so it was giving zero values != possible), implemented new equations...
 *      160628    S.D. Petrovic     Rewrote all transformation angle equations!!
 *      160705    S.D. Petrovic     Added corrections for rounding errors
 *      160712    S.D. Petrovic     Tested new way of writing recurrence relations using Bergsma's thesis
 *
 *    References
 *
 *    Notes
 *
 */





#ifndef AUXILIARY_H
#define AUXILIARY_H


#include <iostream>
#include <iomanip>
#include <cmath>

#include <Eigen/Core>

#include <tudat/Tudat/Mathematics/BasicMathematics/linearAlgebra.h>
#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>

#include <tudatApplications/thesisProject/physicalConstantsUpdated.h>

using namespace std;

class AuxiliaryCartesian

        /* This class will describe the different auxiliary equations, derivatives and functions
         * These are represented by:
         *
         *  - xn                    Auxiliary equation with n is 1, ..., 34
         *  - un                    Auxiliary derivatives n is 1, ..., 34
         *  - wn,m                  Auxiliary functions n is 1, ..., 34 and m is 1, ..., 24
         *
         */





{
public:

    /* In this case, the constructor only takes celestial body and vehicle constant input. The class function will contain the variable parameters
     *
     *    // The diferent celestial body constant parameters and polynomial coefficient parameter matrices
     *
     * const double adiabeticIndex                                   // gamma_a      adiabetic index
     * const double specificGasConstant                        // Rstar    [m^2/(s^2*K)]    specific gas constant
     * const double standardGravitationalParameter  // mu_M     [m^3/s^2]    standard gravitational parameter
     * const double rotationalVelocity                          // rotational velocity of Mars  [rad/s]
     * const double primeMeridianAngle                           // OmegaP   [rad]   relative angle between the prime meridian and the x-axis
     * const double inertialFrameTime                             // t0       [s]    time between the start time and the time that the inertial frame was set
     * const Eigen::MatrixXd temperaturePolyCoefficients // PTn    temperature polynomial coefficients
     * const Eigen::MatrixXd temperatureAltitudeRanges     // altitude range per section for the temperature-altitude curve [km MOLA]
     * const Eigen::VectorXd densityPolyCoefficients         // Prho n density polynomial coefficients
     *
     *     // The differnt vehicle constant parameters and polynomial coefficients
     *
     * const double Thrust                                                         // T   [N]  engine nominal thrust
     * const double specificImpulse                                         // Isp [s]    engine nominal specific impulse
     * const double referenceArea                                            // S [m^2]  vehicle reference area
     * const Eigen::MatrixXd dragCoefficientPolyCoefficients            // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
     * const Eigen::MatrixXd dragCoefficientMachRanges                        // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient
     *
     *
     *
     */


    AuxiliaryCartesian(const double adiabeticIndex_, const double specificGasConstant_, const double standardGravitationalParameter_, const double rotationalVelocity_, const double primeMeridianAngle_,
              const double inertialFrameTime_, const double bodyReferenceRadius_, const Eigen::MatrixXd temperaturePolyCoefficients_, const Eigen::MatrixXd temperatureAltitudeRanges_,
              const Eigen::VectorXd densityPolyCoefficients_, const double Thrust_, const Eigen::MatrixXd thrustAzimuthMatrix_, const Eigen::MatrixXd thrustElevationMatrix_, const double specificImpulse_,
              const double referenceArea_, const Eigen::MatrixXd dragCoefficientPolyCoefficients_, const Eigen::MatrixXd dragCoefficientMachRanges_){

            // Set the different celestial body constant parameters and polynomial coefficient parameter matrices

         adiabeticIndex = adiabeticIndex_;                                   // gamma_a      adiabetic index
         specificGasConstant = specificGasConstant_;                        // Rstar    [m^2/(s^2*K)]    specific gas constant
         standardGravitationalParameter = standardGravitationalParameter_;  // mu_M     [m^3/s^2]    standard gravitational parameter
         rotationalVelocity = rotationalVelocity_;                         // rotational velocity of Mars  [rad/s]
         primeMeridianAngle = primeMeridianAngle_;                          // OmegaP   [rad]   relative angle between the prime meridian and the x-axis
         inertialFrameTime = inertialFrameTime_;                            // t0       [s]    time between the start time and the time that the inertial frame was set
         bodyReferenceRadius = bodyReferenceRadius_;                                  // Rm       [m]     MOLA radius of Mars

         temperaturePolyCoefficients = temperaturePolyCoefficients_; // PTn    temperature polynomial coefficients
         temperatureAltitudeRanges = temperatureAltitudeRanges_;    // altitude range per section for the temperature-altitude curve [km MOLA]
         densityPolyCoefficients = densityPolyCoefficients_;         // Prho n density polynomial coefficients

             // Set the different vehicle constant parameters and polynomial coefficients

         Thrust = Thrust_;                                                         // T   [N]  engine nominal thrust
         thrustAzimuthMatrix = thrustAzimuthMatrix_;                                // psi_T [rad] thrust azimuth angles
         thrustElevationMatrix = thrustElevationMatrix_;                            // epsilon_T [rad] thrust elevation angles
         specificImpulse = specificImpulse_;                                        // Isp [s]    engine nominal specific impulse
         referenceArea = referenceArea_;                                           // S [m^2]  vehicle reference area
         dragCoefficientPolyCoefficients = dragCoefficientPolyCoefficients_;           // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
         dragCoefficientMachRanges = dragCoefficientMachRanges_;                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient


         // Set the flight-path angle and the heading angle to the default values
         FlightPathAngle = 1000;
         HeadingAngle = 1000;           // This is simply done to detect whether or not the value has been set. If the angle does not equal 1000 then the value has been set and has to be used as initial condition

//std::cout<<"verticalInertialFlightPathAngleSet original = "<<verticalInertialFlightPathAngleSet<<std::endl;
//std::cout<<"verticalRotationalFlightPathAngleSet original = "<<verticalRotationalFlightPathAngleSet<<std::endl;



        // Set the booleans to false in case of faulty memory assignment and mistakes in the deletion of the previous class
//       rotationalFlightPathAngleSet = false;         // All of these are used to let the program know that a predefined angle was set and that that angle should be used (initially)
//       inertialFlightPathAngleSet = false;
//       rotationalHeadingAngleSet = false;
//       inertialHeadingAngleSet = false;


//       verticalRotationalFlightPathAngleSet = false;       // All of these are used for the vertical ascent case
//       verticalInertialFlightPathAngleSet = false;
//       verticalRotationalHeadingAngleSet = false;
//       verticalInertialHeadingAngleSet = false;





//std::cout<<"verticalInertialFlightPathAngleSet original 2 = "<<verticalInertialFlightPathAngleSet<<std::endl;




    }

/// Set functions ///

        void setFlightPathAngleAndHeadingAngle(const double FlightPathAngle_ = 1000, const double HeadingAngle_ = 1000){
            if (FlightPathAngle_ >= -tudat::mathematical_constants::LONG_PI/2.0 && FlightPathAngle_ <= tudat::mathematical_constants::LONG_PI/2.0){
                if (HeadingAngle_ >= 0.0 && HeadingAngle_ <= 2.0*tudat::mathematical_constants::LONG_PI){
                    FlightPathAngle = FlightPathAngle_; // Flight path angle in rad
                    HeadingAngle = HeadingAngle_;  // Heading angle in rad

//                    std::cout<<"Flight-path angle and heading angle have been set"<<std::endl;
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






    /// Test functions to revert changes if tolerance is reached



//////////////////////////////////////////////// Auxiliary Equations //////////////////////////////////////////////////////////////////////

    /// Compute the Auxiliary Equations ///

    /* In c++ the convention is that the first entry of a vector is 0, however in this case this first entry will be filled up by w4,2 such that the entry corresponds to the associated auxiliary equation number.
     *
     * The getAuxiliaryEquations function takes three inputs:
     *
     *  - const tudat::basic_mathematics::Vector7d& aState      Which is the current state
     *  - const double time                                     Which is the current time
     *  - const Eigen::Vector3d& thrustAccelerationsBframe      Which are the thrust accelerations in the Bframe
     *
     */

    Eigen::VectorXd getCartesianAuxiliaryEquations( const tudat::basic_mathematics::Vector7d& aState, const double time, const Eigen::Vector3d& thrustAccelerationsBframe){
//        std::cout<<"verticalInertialFlightPathAngleSet eq 1 = "<<verticalInertialFlightPathAngleSet<<std::endl;
//        std::cout<<"verticalRotationalFlightPathAngleSet eq 1 = "<<verticalRotationalFlightPathAngleSet<<std::endl;

        auxiliaryEquationsVector = Eigen::VectorXd::Zero(10);       // Setting the complete vector and filling it with zeros for now



        // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry

        auxiliaryEquationsVector(1) = aState(0);              // x1
        auxiliaryEquationsVector(2) = aState(1);              // x2
        auxiliaryEquationsVector(3) = aState(2);              // x3
        auxiliaryEquationsVector(4) = aState(3);              // x4
        auxiliaryEquationsVector(5) = aState(4);              // x5
        auxiliaryEquationsVector(6) = aState(5);              // x6
        auxiliaryEquationsVector(7) = aState(6);              // x7

       //std::cout<<"Surely this works 1..."<<std::endl;

        auxiliaryEquationsVector(8) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2)+
                auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3) ;              // x8

        auxiliaryEquationsVector(9) = pow(auxiliaryEquationsVector(8), 1.5);              // x9




//        // Please note that the altitude h (x31) is expressed in km MOLA (which is also the input for the density and temperature curves!)
////        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(20)-bodyReferenceRadius)/1000;              // x31 [km]!!!
////        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(20)*1e6-bodyReferenceRadius*1e6)/1e6;              // x31 [km]!!!
//        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(16)-bodyReferenceRadius);              // x31 [km]!!!


//        // Computing the polynomial fit using the altitude and fit parameters for density
//        for (int i = 0; i < 10+1;i++) {

//        auxiliaryEquationsVector(30) += pow(auxiliaryEquationsVector(31),i)*densityPolyCoefficients(i);              // x30
//};

//        // Determine which section of the temperature curve needs to be used and what the corresponding order is
//        // Also, because a computer is less than perfect, a small correction is made to the lower bound of the first section to make sure that the initial altitude is still valid

//        if ((temperatureAltitudeRanges(0,0)-0.000000000001) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(0,1)){

//        sectionT = 0;
//        powerT = 1;

//        }
//        else if (temperatureAltitudeRanges(1,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(1,1)){

//        sectionT = 1;
//        powerT = 3;

//        }
//        else if (temperatureAltitudeRanges(2,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(2,1)){

//        sectionT = 2;
//        powerT = 6;

//        }
//        else if (temperatureAltitudeRanges(3,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(3,1)){

//            sectionT = 3;
//            powerT = 8;
//        }
//        else if (temperatureAltitudeRanges(4,0) <= auxiliaryEquationsVector(31)){

//            sectionT = 4;
//            powerT = 0;
//        }
//        else {


//            std::cerr<<"The current altitude: "<<auxiliaryEquationsVector(31)<<" [km MOLA] is not a valid altitude (lower than the lowest reference altitude)"<<std::endl;

//                       sectionT = 0;
//                        powerT = 1;

//        };
//        //std::cout<<"Surely this works 4..."<<std::endl;

//        // Computing the polynomial fit using the altitude and fit parameters for temperature
//        for (int i=0; i < powerT+1;i++){

//        auxiliaryEquationsVector(34) += pow(auxiliaryEquationsVector(31),i)*temperaturePolyCoefficients(sectionT,i);              // x34

////        std::cout<<"x34 interval "<<i<<" = "<<auxiliaryEquationsVector(34)<<std::endl;

//};

//        auxiliaryEquationsVector(28) = exp(auxiliaryEquationsVector(30));              // x28

//        auxiliaryEquationsVector(33) = sqrt(adiabeticIndex*specificGasConstant*auxiliaryEquationsVector(34));              // x33



//        auxiliaryEquationsVector(32) = auxiliaryEquationsVector(15)/auxiliaryEquationsVector(33);              // x32


//        // Determine which section of the drag coefficient curve needs to be used

//        for (int i=0; i < 5+1; i++){

//            if (dragCoefficientMachRanges(i,0) <= auxiliaryEquationsVector(32) && auxiliaryEquationsVector(32) < dragCoefficientMachRanges(i,1)){

//                sectionCD = i;


//            }


//        };



//        auxiliaryEquationsVector(29) = dragCoefficientPolyCoefficients(sectionCD,1)*auxiliaryEquationsVector(32)+dragCoefficientPolyCoefficients(sectionCD,0);              // x29



//        auxiliaryEquationsVector(27) = 0.5*referenceArea*auxiliaryEquationsVector(28)*auxiliaryEquationsVector(15)*auxiliaryEquationsVector(15)*auxiliaryEquationsVector(29);              // x27


//        auxiliaryEquationsVector(0) = thrustAccelerationsBframe(0)-(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7));              // w4,2





// auxiliaryEquationsVector() = ;               // x


//std::cout<<"Surely this works 2..."<<std::endl;


       return auxiliaryEquationsVector;
    }

//////////////////////////////////////////////// Auxiliary Derivatives //////////////////////////////////////////////////////////////////////

    Eigen::VectorXd getCartesianAuxiliaryDerivatives( const tudat::basic_mathematics::Vector7d& aState, const double time, const Eigen::Vector3d& thrustAccelerationsBframe, const Eigen::VectorXd& auxiliaryEquationsVector){

            //std::cout<<"Surely this works 3..."<<std::endl;

    auxiliaryDerivativesVector = Eigen::VectorXd::Zero(10);       // Setting the complete vector and filling it with zeros for now

    Eigen::MatrixXd auxiliaryFunctionsMatrixDummy = Eigen::MatrixXd::Zero(28,53);       // Setting the complete matrix and filling it with zeros for now

    // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry
    // Which in this case means that the first entry of the vector is 0 and is not used.


//    auxiliaryDerivativesVector(10) = rotationalVelocity;                // u10


    auxiliaryDerivativesVector(1) = auxiliaryEquationsVector(4);                // u1

    auxiliaryDerivativesVector(2) = auxiliaryEquationsVector(5);                // u2

    auxiliaryDerivativesVector(3) = auxiliaryEquationsVector(6);                // u3



    // u4

    auxiliaryFunctionsMatrixDummy(4,1) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2);
    auxiliaryFunctionsMatrixDummy(4,2) = auxiliaryFunctionsMatrixDummy(4,1)+auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3);
    auxiliaryFunctionsMatrixDummy(4,3) = sqrt(auxiliaryFunctionsMatrixDummy(4,2));        // Radius
    auxiliaryFunctionsMatrixDummy(4,4) = sqrt(auxiliaryFunctionsMatrixDummy(4,1));        // 2-D Radius
//    // Avoid singularities
//    if (auxiliaryFunctionsMatrixDummy(4,4) == 0.0){
//        auxiliaryFunctionsMatrixDummy(4,5) = 0.0;
//        auxiliaryFunctionsMatrixDummy(4,6) = 1.0;
//    }
//    else{
//    auxiliaryFunctionsMatrixDummy(4,5) = auxiliaryEquationsVector(2)/auxiliaryFunctionsMatrixDummy(4,4);  // sin(lambda)
//    auxiliaryFunctionsMatrixDummy(4,6) = auxiliaryEquationsVector(1)/auxiliaryFunctionsMatrixDummy(4,4);  // cos(lambda)
//    }
//    auxiliaryFunctionsMatrixDummy(4,7) = auxiliaryEquationsVector(3)/auxiliaryFunctionsMatrixDummy(4,3);  // sin(delta)
//    auxiliaryFunctionsMatrixDummy(4,8) = auxiliaryFunctionsMatrixDummy(4,4)/auxiliaryFunctionsMatrixDummy(4,3);    // cos(delta)
    auxiliaryFunctionsMatrixDummy(4,9) = auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2);
    auxiliaryFunctionsMatrixDummy(4,10) = auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1);
    auxiliaryFunctionsMatrixDummy(4,11) = auxiliaryFunctionsMatrixDummy(4,9)*auxiliaryFunctionsMatrixDummy(4,9)+auxiliaryFunctionsMatrixDummy(4,10)*auxiliaryFunctionsMatrixDummy(4,10)+auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6);
    auxiliaryFunctionsMatrixDummy(4,12)= sqrt(auxiliaryFunctionsMatrixDummy(4,11));


    if (auxiliaryFunctionsMatrixDummy(4,12) == 0.0){
        auxiliaryFunctionsMatrixDummy(4,13) = 0.0;
        auxiliaryFunctionsMatrixDummy(4,14) = 0.0;
        auxiliaryFunctionsMatrixDummy(4,15) = 0.0;
    }
    else{
    auxiliaryFunctionsMatrixDummy(4,13) = auxiliaryFunctionsMatrixDummy(4,9)/auxiliaryFunctionsMatrixDummy(4,12); // X1
    auxiliaryFunctionsMatrixDummy(4,14) = auxiliaryFunctionsMatrixDummy(4,10)/auxiliaryFunctionsMatrixDummy(4,12);  // X2
    auxiliaryFunctionsMatrixDummy(4,15) = auxiliaryEquationsVector(6)/auxiliaryFunctionsMatrixDummy(4,12);  // X3
}

//    auxiliaryFunctionsMatrixDummy(4,16) = auxiliaryFunctionsMatrixDummy(4,10)*auxiliaryEquationsVector(3)-auxiliaryEquationsVector(6)*auxiliaryEquationsVector(2); // Bergsma
//    auxiliaryFunctionsMatrixDummy(4,17) = auxiliaryEquationsVector(6)*auxiliaryEquationsVector(1)-auxiliaryFunctionsMatrixDummy(4,9)*auxiliaryEquationsVector(3);
//    auxiliaryFunctionsMatrixDummy(4,18) = auxiliaryFunctionsMatrixDummy(4,9)*auxiliaryEquationsVector(2)-auxiliaryFunctionsMatrixDummy(4,10)*auxiliaryEquationsVector(1);

    auxiliaryFunctionsMatrixDummy(4,16) = auxiliaryEquationsVector(6)*auxiliaryEquationsVector(2)-auxiliaryFunctionsMatrixDummy(4,10)*auxiliaryEquationsVector(3);
    auxiliaryFunctionsMatrixDummy(4,17) = auxiliaryFunctionsMatrixDummy(4,9)*auxiliaryEquationsVector(3)-auxiliaryEquationsVector(6)*auxiliaryEquationsVector(1);
    auxiliaryFunctionsMatrixDummy(4,18) = auxiliaryFunctionsMatrixDummy(4,10)*auxiliaryEquationsVector(1)-auxiliaryFunctionsMatrixDummy(4,9)*auxiliaryEquationsVector(2);


    auxiliaryFunctionsMatrixDummy(4,19) = auxiliaryFunctionsMatrixDummy(4,16)*auxiliaryFunctionsMatrixDummy(4,16)+auxiliaryFunctionsMatrixDummy(4,17)*auxiliaryFunctionsMatrixDummy(4,17)+auxiliaryFunctionsMatrixDummy(4,18)*auxiliaryFunctionsMatrixDummy(4,18);
    auxiliaryFunctionsMatrixDummy(4,20) = sqrt(auxiliaryFunctionsMatrixDummy(4,19));

    if (auxiliaryFunctionsMatrixDummy(4,12) == 0.0){
        auxiliaryFunctionsMatrixDummy(4,21) = 0.0;
        auxiliaryFunctionsMatrixDummy(4,22) = 0.0;
        auxiliaryFunctionsMatrixDummy(4,23) = 0.0;

    }
    else{
    auxiliaryFunctionsMatrixDummy(4,21) = auxiliaryFunctionsMatrixDummy(4,16)/auxiliaryFunctionsMatrixDummy(4,20); // Y1
    auxiliaryFunctionsMatrixDummy(4,22) = auxiliaryFunctionsMatrixDummy(4,17)/auxiliaryFunctionsMatrixDummy(4,20); // Y2
    auxiliaryFunctionsMatrixDummy(4,23) = auxiliaryFunctionsMatrixDummy(4,18)/auxiliaryFunctionsMatrixDummy(4,20); // Y3
    }
    auxiliaryFunctionsMatrixDummy(4,24) = auxiliaryFunctionsMatrixDummy(4,14)*auxiliaryFunctionsMatrixDummy(4,23)-auxiliaryFunctionsMatrixDummy(4,15)*auxiliaryFunctionsMatrixDummy(4,22); // Z1
    auxiliaryFunctionsMatrixDummy(4,25) = auxiliaryFunctionsMatrixDummy(4,15)*auxiliaryFunctionsMatrixDummy(4,21)-auxiliaryFunctionsMatrixDummy(4,13)*auxiliaryFunctionsMatrixDummy(4,23); // Z2
    auxiliaryFunctionsMatrixDummy(4,40) = auxiliaryFunctionsMatrixDummy(4,13)*auxiliaryFunctionsMatrixDummy(4,22)-auxiliaryFunctionsMatrixDummy(4,14)*auxiliaryFunctionsMatrixDummy(4,21); // Z3


 //std::cout<<"Surely this works 6..."<<std::endl;
//    auxiliaryFunctionsMatrixDummy(4,13) = -auxiliaryFunctionsMatrixDummy(4,6)*auxiliaryFunctionsMatrixDummy(4,7);
//    auxiliaryFunctionsMatrixDummy(4,14) = -auxiliaryFunctionsMatrixDummy(4,7)*auxiliaryFunctionsMatrixDummy(4,5);
//    auxiliaryFunctionsMatrixDummy(4,15) = -auxiliaryFunctionsMatrixDummy(4,8)*auxiliaryFunctionsMatrixDummy(4,6);
//    auxiliaryFunctionsMatrixDummy(4,16) = -auxiliaryFunctionsMatrixDummy(4,8)*auxiliaryFunctionsMatrixDummy(4,5);
//    auxiliaryFunctionsMatrixDummy(4,17) = auxiliaryEquationsVector(6)*auxiliaryFunctionsMatrixDummy(4,8)+auxiliaryFunctionsMatrixDummy(4,9)*auxiliaryFunctionsMatrixDummy(4,13)+auxiliaryFunctionsMatrixDummy(4,10)*auxiliaryFunctionsMatrixDummy(4,14);
//    auxiliaryFunctionsMatrixDummy(4,18) = auxiliaryFunctionsMatrixDummy(4,10)*auxiliaryFunctionsMatrixDummy(4,6)-auxiliaryFunctionsMatrixDummy(4,9)*auxiliaryFunctionsMatrixDummy(4,5);
//    auxiliaryFunctionsMatrixDummy(4,19) = auxiliaryFunctionsMatrixDummy(4,9)*auxiliaryFunctionsMatrixDummy(4,15)-auxiliaryEquationsVector(6)*auxiliaryFunctionsMatrixDummy(4,7)+auxiliaryFunctionsMatrixDummy(4,10)*auxiliaryFunctionsMatrixDummy(4,16);
//    auxiliaryFunctionsMatrixDummy(4,20) = auxiliaryFunctionsMatrixDummy(4,17)*auxiliaryFunctionsMatrixDummy(4,17)+auxiliaryFunctionsMatrixDummy(4,18)*auxiliaryFunctionsMatrixDummy(4,18);
//    auxiliaryFunctionsMatrixDummy(4,21) = sqrt(auxiliaryFunctionsMatrixDummy(4,20));

////    auxiliaryFunctionsMatrixDummy(4,11) = auxiliaryFunctionsMatrixDummy(4,9)*auxiliaryFunctionsMatrixDummy(4,9)+auxiliaryFunctionsMatrixDummy(4,10)*auxiliaryFunctionsMatrixDummy(4,10)+auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6);
////    auxiliaryFunctionsMatrixDummy(4,11) = auxiliaryFunctionsMatrixDummy(4,20)+auxiliaryFunctionsMatrixDummy(4,19)*auxiliaryFunctionsMatrixDummy(4,19);
////    auxiliaryFunctionsMatrixDummy(4,12)= sqrt(auxiliaryFunctionsMatrixDummy(4,11));



//    if (auxiliaryFunctionsMatrixDummy(4,21) == 0.0){
//        auxiliaryFunctionsMatrixDummy(4,22) = sin(HeadingAngle);
//        auxiliaryFunctionsMatrixDummy(4,23) = cos(HeadingAngle);
//        if (abs(auxiliaryFunctionsMatrixDummy(4,22)) < 1.2e-16){
//            auxiliaryFunctionsMatrixDummy(4,22) = 0.0;
//        }
//        if (abs(auxiliaryFunctionsMatrixDummy(4,23)) < 6.2e-17){
//            auxiliaryFunctionsMatrixDummy(4,23) = 0.0;
//        }
//    }
//    else{
//    auxiliaryFunctionsMatrixDummy(4,22) = auxiliaryFunctionsMatrixDummy(4,18)/auxiliaryFunctionsMatrixDummy(4,21);     // sin(chi)
//    auxiliaryFunctionsMatrixDummy(4,23) = auxiliaryFunctionsMatrixDummy(4,17)/auxiliaryFunctionsMatrixDummy(4,21);     // cos(chi)
//}
//     //std::cout<<"Surely this works 7..."<<std::endl;
//    if (auxiliaryFunctionsMatrixDummy(4,12) == 0.0 || time == 0.0){
//        auxiliaryFunctionsMatrixDummy(4,24) = sin(FlightPathAngle);
//        auxiliaryFunctionsMatrixDummy(4,25) = cos(FlightPathAngle);
//        if (abs(auxiliaryFunctionsMatrixDummy(4,24)) < 1.2e-16){
//            auxiliaryFunctionsMatrixDummy(4,24) = 0.0;
//        }
//        if (abs(auxiliaryFunctionsMatrixDummy(4,25)) < 6.2e-17){
//            auxiliaryFunctionsMatrixDummy(4,25) = 0.0;
//        }
//    }
//    else{
//    auxiliaryFunctionsMatrixDummy(4,24) = -auxiliaryFunctionsMatrixDummy(4,19)/auxiliaryFunctionsMatrixDummy(4,12);    // sin(gamma)
//    auxiliaryFunctionsMatrixDummy(4,25) = auxiliaryFunctionsMatrixDummy(4,21)/auxiliaryFunctionsMatrixDummy(4,12);     // cos(gamma)
//  }
    /// Debug ///
//    std::cout<<"uw4,25 = "<<auxiliaryFunctionsMatrixDummy(4,25)<<std::endl;

     //std::cout<<"Surely this works 8..."<<std::endl;
    /// Debug ///



    auxiliaryFunctionsMatrixDummy(27,1) = auxiliaryFunctionsMatrixDummy(4,3) - bodyReferenceRadius; // h
    auxiliaryFunctionsMatrixDummy(27,2) = auxiliaryFunctionsMatrixDummy(27,1)*auxiliaryFunctionsMatrixDummy(27,1);
    auxiliaryFunctionsMatrixDummy(27,3) = pow(auxiliaryFunctionsMatrixDummy(27,1),3);
    auxiliaryFunctionsMatrixDummy(27,4) = pow(auxiliaryFunctionsMatrixDummy(27,1),4);
    auxiliaryFunctionsMatrixDummy(27,5) = pow(auxiliaryFunctionsMatrixDummy(27,1),5);
    auxiliaryFunctionsMatrixDummy(27,6) = pow(auxiliaryFunctionsMatrixDummy(27,1),6);
    auxiliaryFunctionsMatrixDummy(27,7) = pow(auxiliaryFunctionsMatrixDummy(27,1),7);
    auxiliaryFunctionsMatrixDummy(27,8) = pow(auxiliaryFunctionsMatrixDummy(27,1),8);
    auxiliaryFunctionsMatrixDummy(27,9) = pow(auxiliaryFunctionsMatrixDummy(27,1),9);
    auxiliaryFunctionsMatrixDummy(27,10) = pow(auxiliaryFunctionsMatrixDummy(27,1),10);

    // Computing the polynomial fit using the altitude and fit parameters for density
            for (int i = 0; i < 10+1;i++) {

                if (i == 0){
                    auxiliaryFunctionsMatrixDummy(27,11) = densityPolyCoefficients(i);
                }
            else{
            auxiliaryFunctionsMatrixDummy(27,11) += auxiliaryFunctionsMatrixDummy(27,i)*densityPolyCoefficients(i);
    }};

    auxiliaryFunctionsMatrixDummy(27,12) = exp(auxiliaryFunctionsMatrixDummy(27,11)); // Air density

            // Determine which section of the temperature curve needs to be used and what the corresponding order is
            // Also, because a computer is less than perfect, a small correction is made to the lower bound of the first section to make sure that the initial altitude is still valid

            if ((temperatureAltitudeRanges(0,0)-0.000000000001) <= auxiliaryFunctionsMatrixDummy(27,1) && auxiliaryFunctionsMatrixDummy(27,1) < temperatureAltitudeRanges(0,1)){

            sectionT = 0;
            powerT = 1;

            }
            else if (temperatureAltitudeRanges(1,0) <= auxiliaryFunctionsMatrixDummy(27,1) && auxiliaryFunctionsMatrixDummy(27,1) < temperatureAltitudeRanges(1,1)){

            sectionT = 1;
            powerT = 2;

            }
            else if (temperatureAltitudeRanges(2,0) <= auxiliaryFunctionsMatrixDummy(27,1) && auxiliaryFunctionsMatrixDummy(27,1) < temperatureAltitudeRanges(2,1)){

            sectionT = 2;
            powerT = 6;

            }
            else if (temperatureAltitudeRanges(3,0) <= auxiliaryFunctionsMatrixDummy(27,1) && auxiliaryFunctionsMatrixDummy(27,1) < temperatureAltitudeRanges(3,1)){

                sectionT = 3;
                powerT = 8;
            }
            else if (temperatureAltitudeRanges(4,0) <= auxiliaryFunctionsMatrixDummy(27,1)){

                sectionT = 4;
                powerT = 0;
            }
            else {


                std::cerr<<"The current altitude: "<<auxiliaryFunctionsMatrixDummy(27,1)<<" [km MOLA] is not a valid altitude (lower than the lowest reference altitude)"<<std::endl;

                           sectionT = 0;
                            powerT = 1;

            };

                    // Computing the polynomial fit using the altitude and fit parameters for temperature
                    for (int i=0; i < powerT+1;i++){

                        if (i == 0){
                            auxiliaryFunctionsMatrixDummy(27,13) = temperaturePolyCoefficients(sectionT,i);
                        }
                        else {
                    auxiliaryFunctionsMatrixDummy(27,13) += auxiliaryFunctionsMatrixDummy(27,i)*temperaturePolyCoefficients(sectionT,i);              // Air temperature


            }};

                    //std::cout<<"Surely this works 4..."<<std::endl;

    auxiliaryFunctionsMatrixDummy(27,14) = sqrt(adiabeticIndex*specificGasConstant*auxiliaryFunctionsMatrixDummy(27,13)); // Speed of sound
    auxiliaryFunctionsMatrixDummy(27,15) = auxiliaryFunctionsMatrixDummy(4,12)/auxiliaryFunctionsMatrixDummy(27,14); // Mach number

    /// Debug ///
 //std::cout<<"Surely this works 4.5..."<<std::endl;
//  std::cout<<"Temp = "<<auxiliaryFunctionsMatrixDummy(27,13)<<std::endl;
// std::cout<<"a = "<<auxiliaryFunctionsMatrixDummy(27,14)<<std::endl;
// std::cout<<"Mach number = "<<auxiliaryFunctionsMatrixDummy(27,15)<<std::endl;
 /// Debug ///
            // Determine which section of the drag coefficient curve needs to be used

            for (int i=0; i < 5+1; i++){

                if (dragCoefficientMachRanges(i,0) <= auxiliaryFunctionsMatrixDummy(27,15) && auxiliaryFunctionsMatrixDummy(27,15) < dragCoefficientMachRanges(i,1)){

                    sectionCD = i;

//                    std::cout<<"sectionCD = "<<sectionCD<<std::endl;


                }


            };



            auxiliaryFunctionsMatrixDummy(27,16) = dragCoefficientPolyCoefficients(sectionCD,1)*auxiliaryFunctionsMatrixDummy(27,15)+dragCoefficientPolyCoefficients(sectionCD,0);              // Drag coefficient

             //std::cout<<"Surely this works 5..."<<std::endl;

    auxiliaryFunctionsMatrixDummy(27,17) = auxiliaryFunctionsMatrixDummy(4,12)*auxiliaryFunctionsMatrixDummy(4,12);
    auxiliaryFunctionsMatrixDummy(27,18) = auxiliaryFunctionsMatrixDummy(27,17)*auxiliaryFunctionsMatrixDummy(27,16);
    auxiliaryFunctionsMatrixDummy(27,19) = 0.5*referenceArea*auxiliaryFunctionsMatrixDummy(27,18)*auxiliaryFunctionsMatrixDummy(27,12);    // Drag

    // Determining the Thrust azimuth (psiT) and the Thrust elevation (epsilonT) angles based on the altitude
            // Determine the proper azimuth value for the current altitude section
            int sectionThrustAz = 0;    // Set the current azimuth value to the default first section
            for (int i = 0; i < thrustAzimuthMatrix.rows();i++){
                if (thrustAzimuthMatrix(i,0) <= auxiliaryFunctionsMatrixDummy(27,1) && auxiliaryFunctionsMatrixDummy(27,1) < thrustAzimuthMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                    sectionThrustAz = i;
//                        std::cout<<"sectionThrustAz = "<<sectionThrustAz+1<<std::endl;
                }
            }

                const double thrustAzimuth = thrustAzimuthMatrix(sectionThrustAz,2); // Set the thrust azimuth to the current azimuth corresponding to the current altitude section

            // Determine the proper elevation value for the current altitude section
            int sectionThrustEl = 0;    // Set the current elevation value to the default first section
            for (int i = 0; i < thrustElevationMatrix.rows();i++){
                if (thrustElevationMatrix(i,0) <= auxiliaryFunctionsMatrixDummy(27,1) && auxiliaryFunctionsMatrixDummy(27,1) < thrustElevationMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                    sectionThrustEl = i;
                }
            }

                const double thrustElevation = thrustElevationMatrix(sectionThrustEl,2); // Set the thrust elevation to the current elevation corresponding to the current altitude section

//            const double thrustAzimuth = 0.0; // Test
//            const double thrustElevation = 0.0; // Test


            auxiliaryFunctionsMatrixDummy(4,26) = cos(thrustAzimuth);
            auxiliaryFunctionsMatrixDummy(4,27) = cos(thrustElevation);
            auxiliaryFunctionsMatrixDummy(4,28) = sin(thrustAzimuth);
            auxiliaryFunctionsMatrixDummy(4,29) = sin(thrustElevation);
            auxiliaryFunctionsMatrixDummy(4,30) = auxiliaryFunctionsMatrixDummy(4,26)*auxiliaryFunctionsMatrixDummy(4,27);
            auxiliaryFunctionsMatrixDummy(4,31) = auxiliaryFunctionsMatrixDummy(4,27)*auxiliaryFunctionsMatrixDummy(4,28);
            auxiliaryFunctionsMatrixDummy(4,32) = 1/auxiliaryEquationsVector(7);
            auxiliaryFunctionsMatrixDummy(4,33) = Thrust*auxiliaryFunctionsMatrixDummy(4,32);
            auxiliaryFunctionsMatrixDummy(4,34) = auxiliaryFunctionsMatrixDummy(4,33)*auxiliaryFunctionsMatrixDummy(4,30);

    auxiliaryFunctionsMatrixDummy(4,35) = auxiliaryFunctionsMatrixDummy(27,19)/auxiliaryEquationsVector(7);
    auxiliaryFunctionsMatrixDummy(4,36) = auxiliaryFunctionsMatrixDummy(4,34)-auxiliaryFunctionsMatrixDummy(4,35); //
//    auxiliaryFunctionsMatrixDummy(4,36) = 0.0; // Test
    auxiliaryFunctionsMatrixDummy(4,37) = auxiliaryFunctionsMatrixDummy(4,33)*auxiliaryFunctionsMatrixDummy(4,31); //
//    auxiliaryFunctionsMatrixDummy(4,37) = 0.0; // Test
    auxiliaryFunctionsMatrixDummy(4,38) = auxiliaryFunctionsMatrixDummy(4,33)*auxiliaryFunctionsMatrixDummy(4,29); //
//    auxiliaryFunctionsMatrixDummy(4,38) = 0.0; // Test

    auxiliaryFunctionsMatrixDummy(4,39) = -standardGravitationalParameter*auxiliaryEquationsVector(1)/auxiliaryEquationsVector(9);    //
//    auxiliaryFunctionsMatrixDummy(4,40) = -auxiliaryFunctionsMatrixDummy(4,7)*auxiliaryFunctionsMatrixDummy(4,23);
//    auxiliaryFunctionsMatrixDummy(4,41) = auxiliaryFunctionsMatrixDummy(4,8)*auxiliaryFunctionsMatrixDummy(4,24);
//    auxiliaryFunctionsMatrixDummy(4,42) = -auxiliaryFunctionsMatrixDummy(4,5)*auxiliaryFunctionsMatrixDummy(4,22);
//    auxiliaryFunctionsMatrixDummy(4,43) = -auxiliaryFunctionsMatrixDummy(4,5)*auxiliaryFunctionsMatrixDummy(4,23);
//    auxiliaryFunctionsMatrixDummy(4,44) = -auxiliaryFunctionsMatrixDummy(4,8)*auxiliaryFunctionsMatrixDummy(4,25);
//    auxiliaryFunctionsMatrixDummy(4,45) = auxiliaryFunctionsMatrixDummy(4,40)*auxiliaryFunctionsMatrixDummy(4,25);
//    auxiliaryFunctionsMatrixDummy(4,46) = auxiliaryFunctionsMatrixDummy(4,42)*auxiliaryFunctionsMatrixDummy(4,25);
//    auxiliaryFunctionsMatrixDummy(4,47) = -auxiliaryFunctionsMatrixDummy(4,13)*auxiliaryFunctionsMatrixDummy(4,22);
//    auxiliaryFunctionsMatrixDummy(4,48) = auxiliaryFunctionsMatrixDummy(4,40)*auxiliaryFunctionsMatrixDummy(4,24);
//    auxiliaryFunctionsMatrixDummy(4,49) = auxiliaryFunctionsMatrixDummy(4,42)*auxiliaryFunctionsMatrixDummy(4,24);
//    auxiliaryFunctionsMatrixDummy(4,50) = auxiliaryFunctionsMatrixDummy(4,6)*(auxiliaryFunctionsMatrixDummy(4,45)+auxiliaryFunctionsMatrixDummy(4,41))+auxiliaryFunctionsMatrixDummy(4,46);
//    auxiliaryFunctionsMatrixDummy(4,51) = auxiliaryFunctionsMatrixDummy(4,6)*(auxiliaryFunctionsMatrixDummy(4,48)+auxiliaryFunctionsMatrixDummy(4,44))+auxiliaryFunctionsMatrixDummy(4,49);
//    auxiliaryFunctionsMatrixDummy(4,52) = auxiliaryFunctionsMatrixDummy(4,39)+auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,50)+auxiliaryFunctionsMatrixDummy(4,37)*(auxiliaryFunctionsMatrixDummy(4,47)+auxiliaryFunctionsMatrixDummy(4,43))-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,51);

//    auxiliaryDerivativesVector(4) = auxiliaryFunctionsMatrixDummy(4,52); // u4
    auxiliaryDerivativesVector(4) = auxiliaryFunctionsMatrixDummy(4,39)+auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,13)+auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(4,21)-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,24); // u4

//    std::cout<<"u4 = "<<auxiliaryDerivativesVector(4)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,33) = "<<auxiliaryFunctionsMatrixDummy(4,33)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,31) = "<<auxiliaryFunctionsMatrixDummy(4,31)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,27) = "<<auxiliaryFunctionsMatrixDummy(4,27)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,28) = "<<auxiliaryFunctionsMatrixDummy(4,28)<<std::endl;


//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,39) = "<<auxiliaryFunctionsMatrixDummy(4,39)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,13) = "<<auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,13)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(4,21) = "<<auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(4,21)<<std::endl;
//    std::cout<<"-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,24) = "<<-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,24)<<std::endl;

//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,37) = "<<auxiliaryFunctionsMatrixDummy(4,37)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,21) = "<<auxiliaryFunctionsMatrixDummy(4,21)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,32) = "<<auxiliaryFunctionsMatrixDummy(4,32)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,31) = "<<auxiliaryFunctionsMatrixDummy(4,31)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,33) = "<<auxiliaryFunctionsMatrixDummy(4,33)<<std::endl;





    //std::cout<<"Surely this works 9..."<<std::endl;

//    auxiliaryDerivativesVector(4) = -standardGravitationalParameter*(auxiliaryEquationsVector(1)/auxiliaryEquationsVector(9))+auxiliaryEquationsVector(0)*
//            (cx10x11*(-sx12*cx13*cx14+cx12*sx14)-
//                 sx10x11*sx13*cx14)+
//            thrustAccelerationsBframe(1)*(cx10x11*sx12*sx13-sx10x11*cx13)+
//            thrustAccelerationsBframe(2)*(cx10x11*(-sx12*cx13*sx14-
//                                                                                                          cx12*cx14)-
//                                          sx10x11*sx13*sx14);                // u4

    // u5

    auxiliaryFunctionsMatrixDummy(5,1) = -standardGravitationalParameter*auxiliaryEquationsVector(2)/auxiliaryEquationsVector(9);
//    auxiliaryFunctionsMatrixDummy(5,2) = auxiliaryFunctionsMatrixDummy(4,6)*auxiliaryFunctionsMatrixDummy(4,22);
//    auxiliaryFunctionsMatrixDummy(5,3) = auxiliaryFunctionsMatrixDummy(4,5)*(auxiliaryFunctionsMatrixDummy(4,45)+auxiliaryFunctionsMatrixDummy(4,41)+auxiliaryFunctionsMatrixDummy(5,2)*auxiliaryFunctionsMatrixDummy(4,25));
//    auxiliaryFunctionsMatrixDummy(5,4) = -auxiliaryFunctionsMatrixDummy(4,14)*auxiliaryFunctionsMatrixDummy(4,22)+auxiliaryFunctionsMatrixDummy(4,6)*auxiliaryFunctionsMatrixDummy(4,23);
//    auxiliaryFunctionsMatrixDummy(5,5) = auxiliaryFunctionsMatrixDummy(4,5)*(auxiliaryFunctionsMatrixDummy(4,48)+auxiliaryFunctionsMatrixDummy(4,44))+auxiliaryFunctionsMatrixDummy(5,2)*auxiliaryFunctionsMatrixDummy(4,24);
//    auxiliaryFunctionsMatrixDummy(5,6) = auxiliaryFunctionsMatrixDummy(5,1)+auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(5,3)+auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(5,4)-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(5,5);

//    auxiliaryDerivativesVector(5) = auxiliaryFunctionsMatrixDummy(5,6);  // u5
    auxiliaryDerivativesVector(5) = auxiliaryFunctionsMatrixDummy(5,1)+auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,14)+auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(4,22)-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,25); // u5

//    std::cout<<"u5 = "<<auxiliaryDerivativesVector(5)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(5,1) = "<<auxiliaryFunctionsMatrixDummy(5,1)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,14) = "<<auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,14)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(4,22) = "<<auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(4,22)<<std::endl;
//    std::cout<<"-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,25) = "<<-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,25)<<std::endl;



//    auxiliaryDerivativesVector(5) = -standardGravitationalParameter*(auxiliaryEquationsVector(2)/auxiliaryEquationsVector(9))+auxiliaryEquationsVector(0)*
//            (sx10x11*(-sx12*cx13*cx14+cx12*sx14)+
//                 cx10x11*sx13*cx14)+
//            thrustAccelerationsBframe(1)*(sx10x11*sx12*sx13+cx10x11*cx13)+
//            thrustAccelerationsBframe(2)*(sx10x11*(-sx12*cx13*sx14-
//                                                                                                          cx12*cx14)+
//                                          cx10x11*sx13*sx14);                // u5

    // u6

    auxiliaryFunctionsMatrixDummy(6,1) = -standardGravitationalParameter*auxiliaryEquationsVector(3)/auxiliaryEquationsVector(9);
//    auxiliaryFunctionsMatrixDummy(6,2) = auxiliaryFunctionsMatrixDummy(4,7)*auxiliaryFunctionsMatrixDummy(4,24);
//    auxiliaryFunctionsMatrixDummy(6,3) = auxiliaryFunctionsMatrixDummy(4,8)*auxiliaryFunctionsMatrixDummy(4,22);
//    auxiliaryFunctionsMatrixDummy(6,4) = -auxiliaryFunctionsMatrixDummy(4,7)*auxiliaryFunctionsMatrixDummy(4,25);
//    auxiliaryFunctionsMatrixDummy(6,5) = -auxiliaryFunctionsMatrixDummy(4,44)*auxiliaryFunctionsMatrixDummy(4,23)+auxiliaryFunctionsMatrixDummy(6,2);
//    auxiliaryFunctionsMatrixDummy(6,6) = auxiliaryFunctionsMatrixDummy(4,41)*auxiliaryFunctionsMatrixDummy(4,23)+auxiliaryFunctionsMatrixDummy(6,4);
//    auxiliaryFunctionsMatrixDummy(6,7) = auxiliaryFunctionsMatrixDummy(6,1)+auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(6,5)-auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(6,3)-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(6,6);

//    auxiliaryDerivativesVector(6) = auxiliaryFunctionsMatrixDummy(6,7);
    auxiliaryDerivativesVector(6) = auxiliaryFunctionsMatrixDummy(6,1)+auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,15)+auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(4,23)-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,40); // u6


//    std::cout<<"u6 = "<<auxiliaryDerivativesVector(6)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(6,1) = "<<auxiliaryFunctionsMatrixDummy(6,1)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,15) = "<<auxiliaryFunctionsMatrixDummy(4,36)*auxiliaryFunctionsMatrixDummy(4,15)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(4,23) = "<<auxiliaryFunctionsMatrixDummy(4,37)*auxiliaryFunctionsMatrixDummy(4,23)<<std::endl;
//    std::cout<<"-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,40) = "<<-auxiliaryFunctionsMatrixDummy(4,38)*auxiliaryFunctionsMatrixDummy(4,40)<<std::endl;

//    auxiliaryDerivativesVector(6) = -standardGravitationalParameter*(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(9))+auxiliaryEquationsVector(0)*
//            (cx12*cx13*cx14+sx12*sx14)-
//            thrustAccelerationsBframe(1)*cx12*sx13+
//            thrustAccelerationsBframe(2)*(cx12*cx13*sx14-sx12*cx14);                // u6



    auxiliaryDerivativesVector(7) = -(Thrust/(tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION*specificImpulse));                // u7

    auxiliaryDerivativesVector(8) = 2.0*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)+auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6));                // u8

    auxiliaryDerivativesVector(9) = 1.5*((auxiliaryEquationsVector(9)*auxiliaryDerivativesVector(8))/auxiliaryEquationsVector(8));                // u9



// auxiliaryDerivativesVector() = ;                // u


//std::cout<<"Surely this works 10..."<<std::endl;
//std::cout<<"auxiliaryFunctionsMatrixDummy = "<<auxiliaryFunctionsMatrixDummy<<std::endl;
//std::cout<<"auxFunct row 4 = "<<auxiliaryFunctionsMatrixDummy.row(4)<<std::endl;
//std::cout<<"auxFunct row 5 = "<<auxiliaryFunctionsMatrixDummy.row(5)<<std::endl;
//std::cout<<"auxFunct row 6 = "<<auxiliaryFunctionsMatrixDummy.row(6)<<std::endl;
//std::cout<<"auxFunct row 7 = "<<auxiliaryFunctionsMatrixDummy.row(7)<<std::endl;
//std::cout<<"auxFunct row 8 = "<<auxiliaryFunctionsMatrixDummy.row(8)<<std::endl;
//std::cout<<"auxFunct row 9 = "<<auxiliaryFunctionsMatrixDummy.row(9)<<std::endl;
//std::cout<<"auxFunct row 27 = "<<auxiliaryFunctionsMatrixDummy.row(27)<<std::endl;

//std::cout<<"auxiliaryDerivativesVector = "<<auxiliaryDerivativesVector<<std::endl;

//sectionCD = 100;

    return auxiliaryDerivativesVector;
}

//////////////////////////////////////////////// Auxiliary Functions //////////////////////////////////////////////////////////////////////

Eigen::MatrixXd getCartesianAuxiliaryFunctions( const tudat::basic_mathematics::Vector7d& aState, const double time, const Eigen::Vector3d& thrustAccelerationsBframe, const Eigen::VectorXd& auxiliaryEquationsVector,
                                         const Eigen::VectorXd& auxiliaryDerivativesVectorInput){

    auxiliaryFunctionsMatrix = Eigen::MatrixXd::Zero(28,53);       // Setting the complete matrix and filling it with zeros for now

//    Eigen::MatrixXd thrustAzimuthMatrix = MAV.thrustAzimuth();      // Setting the thrust angle matrices
//    Eigen::MatrixXd thrustElevationMatrix = MAV.thrustElevation();

    // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry
    // Which in this case means that the not all positions in the matrix will be used. The other values will simply be 0.

    //std::cout<<"I mean, come on... it works 1"<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrix beginning = "<<auxiliaryFunctionsMatrix<<std::endl;

    // u4

    auxiliaryFunctionsMatrix(4,1) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2);
    auxiliaryFunctionsMatrix(4,2) = auxiliaryFunctionsMatrix(4,1)+auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3);
    auxiliaryFunctionsMatrix(4,3) = sqrt(auxiliaryFunctionsMatrix(4,2));        // Radius
    auxiliaryFunctionsMatrix(4,4) = sqrt(auxiliaryFunctionsMatrix(4,1));        // 2-D Radius

    // Avoid singularities
    if (auxiliaryFunctionsMatrix(4,4) == 0.0){
        auxiliaryFunctionsMatrix(4,5) = 0.0;
        auxiliaryFunctionsMatrix(4,6) = 1.0;
    }
    else{
    auxiliaryFunctionsMatrix(4,5) = auxiliaryEquationsVector(2)/auxiliaryFunctionsMatrix(4,4);  // sin(lambda)
    auxiliaryFunctionsMatrix(4,6) = auxiliaryEquationsVector(1)/auxiliaryFunctionsMatrix(4,4);  // cos(lambda)
   }
    auxiliaryFunctionsMatrix(4,7) = auxiliaryEquationsVector(3)/auxiliaryFunctionsMatrix(4,3);  // sin(delta)
    auxiliaryFunctionsMatrix(4,8) = auxiliaryFunctionsMatrix(4,4)/auxiliaryFunctionsMatrix(4,3);    // cos(delta)
    auxiliaryFunctionsMatrix(4,9) = auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2);
    auxiliaryFunctionsMatrix(4,10) = auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1);

        auxiliaryFunctionsMatrix(4,11) = auxiliaryFunctionsMatrix(4,9)*auxiliaryFunctionsMatrix(4,9)+auxiliaryFunctionsMatrix(4,10)*auxiliaryFunctionsMatrix(4,10)+auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6);
//        auxiliaryFunctionsMatrix(4,11) = auxiliaryFunctionsMatrix(4,20)+auxiliaryFunctionsMatrix(4,19)*auxiliaryFunctionsMatrix(4,19);
        auxiliaryFunctionsMatrix(4,12)= sqrt(auxiliaryFunctionsMatrix(4,11)); // Ground Velocity

        if (auxiliaryFunctionsMatrix(4,12) == 0.0){
            auxiliaryFunctionsMatrix(4,13) = 0.0;
            auxiliaryFunctionsMatrix(4,14) = 0.0;
            auxiliaryFunctionsMatrix(4,15) = 0.0;
        }
        else{

    auxiliaryFunctionsMatrix(4,13) = auxiliaryFunctionsMatrix(4,9)/auxiliaryFunctionsMatrix(4,12); // X1
    auxiliaryFunctionsMatrix(4,14) = auxiliaryFunctionsMatrix(4,10)/auxiliaryFunctionsMatrix(4,12);  // X2
    auxiliaryFunctionsMatrix(4,15) = auxiliaryEquationsVector(6)/auxiliaryFunctionsMatrix(4,12);  // X3

}

//        auxiliaryFunctionsMatrix(4,16) = auxiliaryFunctionsMatrix(4,10)*auxiliaryEquationsVector(3)-auxiliaryEquationsVector(6)*auxiliaryEquationsVector(2); // Bergsma
//        auxiliaryFunctionsMatrix(4,17) = auxiliaryEquationsVector(6)*auxiliaryEquationsVector(1)-auxiliaryFunctionsMatrix(4,9)*auxiliaryEquationsVector(3);
//        auxiliaryFunctionsMatrix(4,18) = auxiliaryFunctionsMatrix(4,9)*auxiliaryEquationsVector(2)-auxiliaryFunctionsMatrix(4,10)*auxiliaryEquationsVector(1);
        auxiliaryFunctionsMatrix(4,16) = auxiliaryEquationsVector(6)*auxiliaryEquationsVector(2)-auxiliaryFunctionsMatrix(4,10)*auxiliaryEquationsVector(3);
        auxiliaryFunctionsMatrix(4,17) = auxiliaryFunctionsMatrix(4,9)*auxiliaryEquationsVector(3)-auxiliaryEquationsVector(6)*auxiliaryEquationsVector(1);
        auxiliaryFunctionsMatrix(4,18) = auxiliaryFunctionsMatrix(4,10)*auxiliaryEquationsVector(1)-auxiliaryFunctionsMatrix(4,9)*auxiliaryEquationsVector(2);


        auxiliaryFunctionsMatrix(4,19) = auxiliaryFunctionsMatrix(4,16)*auxiliaryFunctionsMatrix(4,16)+auxiliaryFunctionsMatrix(4,17)*auxiliaryFunctionsMatrix(4,17)+auxiliaryFunctionsMatrix(4,18)*auxiliaryFunctionsMatrix(4,18);
        auxiliaryFunctionsMatrix(4,20) = sqrt(auxiliaryFunctionsMatrix(4,19));

        if (auxiliaryFunctionsMatrix(4,12) == 0.0){
            auxiliaryFunctionsMatrix(4,21) = 0.0;
            auxiliaryFunctionsMatrix(4,22) = 0.0;
            auxiliaryFunctionsMatrix(4,23) = 0.0;

        }
        else{
        auxiliaryFunctionsMatrix(4,21) = auxiliaryFunctionsMatrix(4,16)/auxiliaryFunctionsMatrix(4,20); // Y1
        auxiliaryFunctionsMatrix(4,22) = auxiliaryFunctionsMatrix(4,17)/auxiliaryFunctionsMatrix(4,20); // Y2
        auxiliaryFunctionsMatrix(4,23) = auxiliaryFunctionsMatrix(4,18)/auxiliaryFunctionsMatrix(4,20); // Y3
        }
        auxiliaryFunctionsMatrix(4,24) = auxiliaryFunctionsMatrix(4,14)*auxiliaryFunctionsMatrix(4,23)-auxiliaryFunctionsMatrix(4,15)*auxiliaryFunctionsMatrix(4,22); // Z1
        auxiliaryFunctionsMatrix(4,25) = auxiliaryFunctionsMatrix(4,15)*auxiliaryFunctionsMatrix(4,21)-auxiliaryFunctionsMatrix(4,13)*auxiliaryFunctionsMatrix(4,23); // Z2
        auxiliaryFunctionsMatrix(4,40) = auxiliaryFunctionsMatrix(4,13)*auxiliaryFunctionsMatrix(4,22)-auxiliaryFunctionsMatrix(4,14)*auxiliaryFunctionsMatrix(4,21); // Z3

    /// Debug ///
//        std::cout<<"auxiliaryFunctionsMatrix(4,16) = "<<auxiliaryFunctionsMatrix(4,16)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,17) = "<<auxiliaryFunctionsMatrix(4,17)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,18) = "<<auxiliaryFunctionsMatrix(4,18)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,19) = "<<auxiliaryFunctionsMatrix(4,19)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,20) = "<<auxiliaryFunctionsMatrix(4,20)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,13) = "<<auxiliaryFunctionsMatrix(4,13)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,14) = "<<auxiliaryFunctionsMatrix(4,14)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,15) = "<<auxiliaryFunctionsMatrix(4,15)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,21) = "<<auxiliaryFunctionsMatrix(4,21)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,22) = "<<auxiliaryFunctionsMatrix(4,22)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,23) = "<<auxiliaryFunctionsMatrix(4,23)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,24) = "<<auxiliaryFunctionsMatrix(4,24)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,25) = "<<auxiliaryFunctionsMatrix(4,25)<<std::endl;
//        std::cout<<"auxiliaryFunctionsMatrix(4,40) = "<<auxiliaryFunctionsMatrix(4,40)<<std::endl;


//    std::cout<<"x1-x2 = "<<auxiliaryEquationsVector(1)-auxiliaryEquationsVector(2)<<std::endl;
//        std::cout<<"w4,4 = "<<auxiliaryFunctionsMatrix(4,4)<<std::endl;
//    std::cout<<"w4,5 = "<<auxiliaryFunctionsMatrix(4,5)<<std::endl;
//    std::cout<<"w4,6 = "<<auxiliaryFunctionsMatrix(4,6)<<std::endl;
//    std::cout<<"w4,7 = "<<auxiliaryFunctionsMatrix(4,7)<<std::endl;
//    std::cout<<"w4,8 = "<<auxiliaryFunctionsMatrix(4,8)<<std::endl;

    /// Debug ///
      //std::cout<<"I mean, come on... it works 2"<<std::endl;




//    auxiliaryFunctionsMatrix(4,13) = -auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,7);
//    auxiliaryFunctionsMatrix(4,14) = -auxiliaryFunctionsMatrix(4,7)*auxiliaryFunctionsMatrix(4,5);
//    auxiliaryFunctionsMatrix(4,15) = -auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,6);
//    auxiliaryFunctionsMatrix(4,16) = -auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,5);
//    auxiliaryFunctionsMatrix(4,17) = auxiliaryEquationsVector(6)*auxiliaryFunctionsMatrix(4,8)+auxiliaryFunctionsMatrix(4,9)*auxiliaryFunctionsMatrix(4,13)+auxiliaryFunctionsMatrix(4,10)*auxiliaryFunctionsMatrix(4,14); // Vx_v
//    auxiliaryFunctionsMatrix(4,18) = auxiliaryFunctionsMatrix(4,10)*auxiliaryFunctionsMatrix(4,6)-auxiliaryFunctionsMatrix(4,9)*auxiliaryFunctionsMatrix(4,5);  // Vy_v
//    auxiliaryFunctionsMatrix(4,19) = auxiliaryFunctionsMatrix(4,9)*auxiliaryFunctionsMatrix(4,15)-auxiliaryEquationsVector(6)*auxiliaryFunctionsMatrix(4,7)+auxiliaryFunctionsMatrix(4,10)*auxiliaryFunctionsMatrix(4,16); // Vz_v
//    auxiliaryFunctionsMatrix(4,20) = auxiliaryFunctionsMatrix(4,17)*auxiliaryFunctionsMatrix(4,17)+auxiliaryFunctionsMatrix(4,18)*auxiliaryFunctionsMatrix(4,18);
//    auxiliaryFunctionsMatrix(4,21) = sqrt(auxiliaryFunctionsMatrix(4,20));

//////    auxiliaryFunctionsMatrix(4,11) = auxiliaryFunctionsMatrix(4,9)*auxiliaryFunctionsMatrix(4,9)+auxiliaryFunctionsMatrix(4,10)*auxiliaryFunctionsMatrix(4,10)+auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6);
////    auxiliaryFunctionsMatrix(4,11) = auxiliaryFunctionsMatrix(4,20)+auxiliaryFunctionsMatrix(4,19)*auxiliaryFunctionsMatrix(4,19);
////    auxiliaryFunctionsMatrix(4,12)= sqrt(auxiliaryFunctionsMatrix(4,11)); // Ground Velocity

////    std::cout<<"w4,19 (Vz_v) = "<<auxiliaryFunctionsMatrix(4,19)<<std::endl;

//    if (auxiliaryFunctionsMatrix(4,21) == 0.0){
//        auxiliaryFunctionsMatrix(4,22) = sin(HeadingAngle);
//        auxiliaryFunctionsMatrix(4,23) = cos(HeadingAngle);
//        if (abs(auxiliaryFunctionsMatrix(4,22)) < 1.2e-16){
//            auxiliaryFunctionsMatrix(4,22) = 0.0;
//        }
//        if (abs(auxiliaryFunctionsMatrix(4,23)) < 6.2e-17){
//            auxiliaryFunctionsMatrix(4,23) = 0.0;
//        }
//    }
//    else{
//    auxiliaryFunctionsMatrix(4,22) = auxiliaryFunctionsMatrix(4,18)/auxiliaryFunctionsMatrix(4,21);     // sin(chi)
//    auxiliaryFunctionsMatrix(4,23) = auxiliaryFunctionsMatrix(4,17)/auxiliaryFunctionsMatrix(4,21);     // cos(chi)
//}

//    if (auxiliaryFunctionsMatrix(4,12) == 0.0 || time == 0.0){
//        auxiliaryFunctionsMatrix(4,24) = sin(FlightPathAngle);
//        auxiliaryFunctionsMatrix(4,25) = cos(FlightPathAngle);
//        if (abs(auxiliaryFunctionsMatrix(4,24)) < 1.2e-16){
//            auxiliaryFunctionsMatrix(4,24) = 0.0;
//        }
//        if (abs(auxiliaryFunctionsMatrix(4,25)) < 6.2e-17){
//            auxiliaryFunctionsMatrix(4,25) = 0.0;
//        }
//    }
//    else{
//    auxiliaryFunctionsMatrix(4,24) = -auxiliaryFunctionsMatrix(4,19)/auxiliaryFunctionsMatrix(4,12);    // sin(gamma)
//    auxiliaryFunctionsMatrix(4,25) = auxiliaryFunctionsMatrix(4,21)/auxiliaryFunctionsMatrix(4,12);     // cos(gamma)
//  }

//    /// Debug ///
////    std::cout<<"w4,25 = "<<auxiliaryFunctionsMatrix(4,25)<<std::endl;
//    std::cout<<"-auxiliaryFunctionsMatrix(4,19) = "<<-auxiliaryFunctionsMatrix(4,19)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrix(4,12) = "<<auxiliaryFunctionsMatrix(4,12)<<std::endl;
//    std::cout<<"tan(gamma) = "<<-auxiliaryFunctionsMatrix(4,19)/auxiliaryFunctionsMatrix(4,21)<<std::endl;
//    std::cout<<"sin(gamma) = tan(gamma)*cos(gamma) = "<<(-auxiliaryFunctionsMatrix(4,19)/auxiliaryFunctionsMatrix(4,21))*auxiliaryFunctionsMatrix(4,25)<<std::endl;
//    /// Debug ///



    auxiliaryFunctionsMatrix(27,1) = auxiliaryFunctionsMatrix(4,3) - bodyReferenceRadius;
    auxiliaryFunctionsMatrix(27,2) = auxiliaryFunctionsMatrix(27,1)*auxiliaryFunctionsMatrix(27,1);
    auxiliaryFunctionsMatrix(27,3) = pow(auxiliaryFunctionsMatrix(27,1),3);
    auxiliaryFunctionsMatrix(27,4) = pow(auxiliaryFunctionsMatrix(27,1),4);
    auxiliaryFunctionsMatrix(27,5) = pow(auxiliaryFunctionsMatrix(27,1),5);
    auxiliaryFunctionsMatrix(27,6) = pow(auxiliaryFunctionsMatrix(27,1),6);
    auxiliaryFunctionsMatrix(27,7) = pow(auxiliaryFunctionsMatrix(27,1),7);
    auxiliaryFunctionsMatrix(27,8) = pow(auxiliaryFunctionsMatrix(27,1),8);
    auxiliaryFunctionsMatrix(27,9) = pow(auxiliaryFunctionsMatrix(27,1),9);
    auxiliaryFunctionsMatrix(27,10) = pow(auxiliaryFunctionsMatrix(27,1),10);

    // Computing the polynomial fit using the altitude and fit parameters for density
            for (int i = 0; i < 10+1;i++) {

                if (i == 0){
                    auxiliaryFunctionsMatrix(27,11) = densityPolyCoefficients(i);
                }
            else{
            auxiliaryFunctionsMatrix(27,11) += auxiliaryFunctionsMatrix(27,i)*densityPolyCoefficients(i);
    }};

    auxiliaryFunctionsMatrix(27,12) = exp(auxiliaryFunctionsMatrix(27,11)); // Air density

            // Determine which section of the temperature curve needs to be used and what the corresponding order is
            // Also, because a computer is less than perfect, a small correction is made to the lower bound of the first section to make sure that the initial altitude is still valid

//            if ((temperatureAltitudeRanges(0,0)-0.000000000001) <= auxiliaryFunctionsMatrix(27,1) && auxiliaryFunctionsMatrix(27,1) < temperatureAltitudeRanges(0,1)){

//            sectionT = 0;
//            powerT = 1;

//            }
//            else if (temperatureAltitudeRanges(1,0) <= auxiliaryFunctionsMatrix(27,1) && auxiliaryFunctionsMatrix(27,1) < temperatureAltitudeRanges(1,1)){

//            sectionT = 1;
//            powerT = 2;

//            }
//            else if (temperatureAltitudeRanges(2,0) <= auxiliaryFunctionsMatrix(27,1) && auxiliaryFunctionsMatrix(27,1) < temperatureAltitudeRanges(2,1)){

//            sectionT = 2;
//            powerT = 6;

//            }
//            else if (temperatureAltitudeRanges(3,0) <= auxiliaryFunctionsMatrix(27,1) && auxiliaryFunctionsMatrix(27,1) < temperatureAltitudeRanges(3,1)){

//                sectionT = 3;
//                powerT = 8;
//            }
//            else if (temperatureAltitudeRanges(4,0) <= auxiliaryFunctionsMatrix(27,1)){

//                sectionT = 4;
//                powerT = 0;
//            }
//            else {


//                std::cerr<<"The current altitude: "<<auxiliaryFunctionsMatrix(27,1)<<" [km MOLA] is not a valid altitude (lower than the lowest reference altitude)"<<std::endl;

//                           sectionT = 0;
//                            powerT = 1;

//            };

                    // Computing the polynomial fit using the altitude and fit parameters for temperature
                    for (int i=0; i < powerT+1;i++){

                        if (i == 0){
                            auxiliaryFunctionsMatrix(27,13) = temperaturePolyCoefficients(sectionT,i);
                        }
                        else {
                    auxiliaryFunctionsMatrix(27,13) += auxiliaryFunctionsMatrix(27,i)*temperaturePolyCoefficients(sectionT,i);              // Air temperature


            }};

    auxiliaryFunctionsMatrix(27,14) = sqrt(adiabeticIndex*specificGasConstant*auxiliaryFunctionsMatrix(27,13)); // Speed of sound
    auxiliaryFunctionsMatrix(27,15) = auxiliaryFunctionsMatrix(4,12)/auxiliaryFunctionsMatrix(27,14); // Mach number

            // Determine which section of the drag coefficient curve needs to be used
//    std::cout<<"sectionCD before = "<<sectionCD<<std::endl;

//            for (int i=0; i < 5+1; i++){

//                if (dragCoefficientMachRanges(i,0) <= auxiliaryFunctionsMatrix(27,15) && auxiliaryFunctionsMatrix(27,15) < dragCoefficientMachRanges(i,1)){

//                    sectionCD = i;

//                    std::cout<<"sectionCD after = "<<sectionCD<<std::endl;


//                }


//            };

  //std::cout<<"I mean, come on... it works 3"<<std::endl;

            auxiliaryFunctionsMatrix(27,16) = dragCoefficientPolyCoefficients(sectionCD,1)*auxiliaryFunctionsMatrix(27,15)+dragCoefficientPolyCoefficients(sectionCD,0);              // Drag coefficient


    auxiliaryFunctionsMatrix(27,17) = auxiliaryFunctionsMatrix(4,12)*auxiliaryFunctionsMatrix(4,12);
    auxiliaryFunctionsMatrix(27,18) = auxiliaryFunctionsMatrix(27,17)*auxiliaryFunctionsMatrix(27,16);
    auxiliaryFunctionsMatrix(27,19) = 0.5*referenceArea*auxiliaryFunctionsMatrix(27,18)*auxiliaryFunctionsMatrix(27,12);    // Drag


    // Determining the Thrust azimuth (psiT) and the Thrust elevation (epsilonT) angles based on the altitude
            // Determine the proper azimuth value for the current altitude section
            int sectionThrustAz = 0;    // Set the current azimuth value to the default first section
            for (int i = 0; i < thrustAzimuthMatrix.rows();i++){
                if (thrustAzimuthMatrix(i,0) <= auxiliaryFunctionsMatrix(27,1) && auxiliaryFunctionsMatrix(27,1) < thrustAzimuthMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                    sectionThrustAz = i;
//                        std::cout<<"sectionThrustAz = "<<sectionThrustAz+1<<std::endl;
                }
            }

                const double thrustAzimuth = thrustAzimuthMatrix(sectionThrustAz,2); // Set the thrust azimuth to the current azimuth corresponding to the current altitude section

            // Determine the proper elevation value for the current altitude section
            int sectionThrustEl = 0;    // Set the current elevation value to the default first section
            for (int i = 0; i < thrustElevationMatrix.rows();i++){
                if (thrustElevationMatrix(i,0) <= auxiliaryFunctionsMatrix(27,1) && auxiliaryFunctionsMatrix(27,1) < thrustElevationMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                    sectionThrustEl = i;
                }
            }

                const double thrustElevation = thrustElevationMatrix(sectionThrustEl,2); // Set the thrust elevation to the current elevation corresponding to the current altitude section

//            const double thrustAzimuth = 0.0; // Test
//            const double thrustElevation = 0.0; // Test

            auxiliaryFunctionsMatrix(4,26) = cos(thrustAzimuth);
            auxiliaryFunctionsMatrix(4,27) = cos(thrustElevation);
            auxiliaryFunctionsMatrix(4,28) = sin(thrustAzimuth);
            auxiliaryFunctionsMatrix(4,29) = sin(thrustElevation);
            auxiliaryFunctionsMatrix(4,30) = auxiliaryFunctionsMatrix(4,26)*auxiliaryFunctionsMatrix(4,27);
            auxiliaryFunctionsMatrix(4,31) = auxiliaryFunctionsMatrix(4,27)*auxiliaryFunctionsMatrix(4,28);
            auxiliaryFunctionsMatrix(4,32) = 1/auxiliaryEquationsVector(7);
            auxiliaryFunctionsMatrix(4,33) = Thrust*auxiliaryFunctionsMatrix(4,32);
            auxiliaryFunctionsMatrix(4,34) = auxiliaryFunctionsMatrix(4,33)*auxiliaryFunctionsMatrix(4,30);
    auxiliaryFunctionsMatrix(4,35) = auxiliaryFunctionsMatrix(27,19)/auxiliaryEquationsVector(7);

    /// Debug ///
//    std::cout<<"w4,35 = "<<auxiliaryFunctionsMatrix(4,35)<<std::endl;
//    std::cout<<"x7 = "<<auxiliaryEquationsVector(7)<<std::endl;
//    std::cout<<"w17,19 = "<<auxiliaryFunctionsMatrix(27,19)<<std::endl;
    /// Debug ///

    auxiliaryFunctionsMatrix(4,36) = auxiliaryFunctionsMatrix(4,34)-auxiliaryFunctionsMatrix(4,35); //
//    auxiliaryFunctionsMatrix(4,36) = 0.0; // Test

    auxiliaryFunctionsMatrix(4,37) = auxiliaryFunctionsMatrix(4,33)*auxiliaryFunctionsMatrix(4,31); //
//    auxiliaryFunctionsMatrix(4,37) = 0.0; // Test

    auxiliaryFunctionsMatrix(4,38) = auxiliaryFunctionsMatrix(4,33)*auxiliaryFunctionsMatrix(4,29); //
//    auxiliaryFunctionsMatrix(4,38) = 0.0; // Test

    auxiliaryFunctionsMatrix(4,39) = -standardGravitationalParameter*auxiliaryEquationsVector(1)/auxiliaryEquationsVector(9);    //
//    auxiliaryFunctionsMatrix(4,40) = -auxiliaryFunctionsMatrix(4,7)*auxiliaryFunctionsMatrix(4,23);
//    auxiliaryFunctionsMatrix(4,41) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,24);
//    auxiliaryFunctionsMatrix(4,42) = -auxiliaryFunctionsMatrix(4,5)*auxiliaryFunctionsMatrix(4,22);
//    auxiliaryFunctionsMatrix(4,43) = -auxiliaryFunctionsMatrix(4,5)*auxiliaryFunctionsMatrix(4,23);
//    auxiliaryFunctionsMatrix(4,44) = -auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,25);
//    auxiliaryFunctionsMatrix(4,45) = auxiliaryFunctionsMatrix(4,40)*auxiliaryFunctionsMatrix(4,25);
//    auxiliaryFunctionsMatrix(4,46) = auxiliaryFunctionsMatrix(4,42)*auxiliaryFunctionsMatrix(4,25);
//    auxiliaryFunctionsMatrix(4,47) = -auxiliaryFunctionsMatrix(4,13)*auxiliaryFunctionsMatrix(4,22);
//    auxiliaryFunctionsMatrix(4,48) = auxiliaryFunctionsMatrix(4,40)*auxiliaryFunctionsMatrix(4,24);
//    auxiliaryFunctionsMatrix(4,49) = auxiliaryFunctionsMatrix(4,42)*auxiliaryFunctionsMatrix(4,24);
//    auxiliaryFunctionsMatrix(4,50) = auxiliaryFunctionsMatrix(4,6)*(auxiliaryFunctionsMatrix(4,45)+auxiliaryFunctionsMatrix(4,41))+auxiliaryFunctionsMatrix(4,46);
//    auxiliaryFunctionsMatrix(4,51) = auxiliaryFunctionsMatrix(4,6)*(auxiliaryFunctionsMatrix(4,48)+auxiliaryFunctionsMatrix(4,44))+auxiliaryFunctionsMatrix(4,49);
//    auxiliaryFunctionsMatrix(4,52) = auxiliaryFunctionsMatrix(4,39)+auxiliaryFunctionsMatrix(4,36)*auxiliaryFunctionsMatrix(4,50)+auxiliaryFunctionsMatrix(4,37)*(auxiliaryFunctionsMatrix(4,47)+auxiliaryFunctionsMatrix(4,43))-auxiliaryFunctionsMatrix(4,38)*auxiliaryFunctionsMatrix(4,51);


    /// Debug ///
//    std::cout<<"auxiliaryFunctionsMatrix(4,52) = "<<auxiliaryFunctionsMatrix(4,52)<<std::endl;
////    std::cout<<"auxiliaryFunctionsMatrix(4,39) = "<<auxiliaryFunctionsMatrix(4,39)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrix(4,36)*auxiliaryFunctionsMatrix(4,50) = "<<auxiliaryFunctionsMatrix(4,36)*auxiliaryFunctionsMatrix(4,50)<<std::endl;
////    std::cout<<"auxiliaryFunctionsMatrix(4,37)*(auxiliaryFunctionsMatrix(4,47)+auxiliaryFunctionsMatrix(4,43)) = "<<auxiliaryFunctionsMatrix(4,37)*(auxiliaryFunctionsMatrix(4,47)+auxiliaryFunctionsMatrix(4,43))<<std::endl;
////    std::cout<<"-auxiliaryFunctionsMatrix(4,38)*auxiliaryFunctionsMatrix(4,51) = "<<-auxiliaryFunctionsMatrix(4,38)*auxiliaryFunctionsMatrix(4,51)<<std::endl;
////    std::cout<<"auxiliaryFunctionsMatrix(4,36) = "<<auxiliaryFunctionsMatrix(4,36)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrix(4,50) = "<<auxiliaryFunctionsMatrix(4,50)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrix(4,6)*(auxiliaryFunctionsMatrix(4,43)+auxiliaryFunctionsMatrix(4,41)) = "<<auxiliaryFunctionsMatrix(4,6)*(auxiliaryFunctionsMatrix(4,43)+auxiliaryFunctionsMatrix(4,41))<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrix(4,46) = "<<auxiliaryFunctionsMatrix(4,46)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrix(4,6) = "<<auxiliaryFunctionsMatrix(4,6)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrix(4,43) = "<<auxiliaryFunctionsMatrix(4,43)<<std::endl;
//    std::cout<<"auxiliaryFunctionsMatrix(4,41) = "<<auxiliaryFunctionsMatrix(4,41)<<std::endl;


    /// Debug ///

    // w5

    auxiliaryFunctionsMatrix(5,1) = -standardGravitationalParameter*auxiliaryEquationsVector(2)/auxiliaryEquationsVector(9);
//    auxiliaryFunctionsMatrix(5,2) = auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,22);
//    auxiliaryFunctionsMatrix(5,3) = auxiliaryFunctionsMatrix(4,5)*(auxiliaryFunctionsMatrix(4,45)+auxiliaryFunctionsMatrix(4,41)+auxiliaryFunctionsMatrix(5,2)*auxiliaryFunctionsMatrix(4,25));
//    auxiliaryFunctionsMatrix(5,4) = -auxiliaryFunctionsMatrix(4,14)*auxiliaryFunctionsMatrix(4,22)+auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,23);
//    auxiliaryFunctionsMatrix(5,5) = auxiliaryFunctionsMatrix(4,5)*(auxiliaryFunctionsMatrix(4,48)+auxiliaryFunctionsMatrix(4,44))+auxiliaryFunctionsMatrix(5,2)*auxiliaryFunctionsMatrix(4,24);
//    auxiliaryFunctionsMatrix(5,6) = auxiliaryFunctionsMatrix(5,1)+auxiliaryFunctionsMatrix(4,36)*auxiliaryFunctionsMatrix(5,3)+auxiliaryFunctionsMatrix(4,37)*auxiliaryFunctionsMatrix(5,4)-auxiliaryFunctionsMatrix(4,38)*auxiliaryFunctionsMatrix(5,5);

  //std::cout<<"I mean, come on... it works 4"<<std::endl;
    // w6

    auxiliaryFunctionsMatrix(6,1) = -standardGravitationalParameter*auxiliaryEquationsVector(3)/auxiliaryEquationsVector(9);
//    auxiliaryFunctionsMatrix(6,2) = auxiliaryFunctionsMatrix(4,7)*auxiliaryFunctionsMatrix(4,24);
//    auxiliaryFunctionsMatrix(6,3) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,22);
//    auxiliaryFunctionsMatrix(6,4) = -auxiliaryFunctionsMatrix(4,7)*auxiliaryFunctionsMatrix(4,25);
//    auxiliaryFunctionsMatrix(6,5) = -auxiliaryFunctionsMatrix(4,44)*auxiliaryFunctionsMatrix(4,23)+auxiliaryFunctionsMatrix(6,2);
//    auxiliaryFunctionsMatrix(6,6) = auxiliaryFunctionsMatrix(4,41)*auxiliaryFunctionsMatrix(4,23)+auxiliaryFunctionsMatrix(6,4);
//    auxiliaryFunctionsMatrix(6,7) = auxiliaryFunctionsMatrix(6,1)+auxiliaryFunctionsMatrix(4,36)*auxiliaryFunctionsMatrix(6,5)-auxiliaryFunctionsMatrix(4,37)*auxiliaryFunctionsMatrix(6,3)-auxiliaryFunctionsMatrix(4,38)*auxiliaryFunctionsMatrix(6,6);



    // w8
    auxiliaryFunctionsMatrix(8,1) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4);
    auxiliaryFunctionsMatrix(8,2) = auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5);
    auxiliaryFunctionsMatrix(8,3) = auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6);





    // w9
    auxiliaryFunctionsMatrix(9,1) = (auxiliaryEquationsVector(9)*auxiliaryDerivativesVectorInput(8))/auxiliaryEquationsVector(8);


// auxiliaryFunctionsMatrix(,)

      //std::cout<<"I mean, come on... it works 5"<<std::endl;
//      std::cout<<"auxiliaryFunctionsMatrix end = "<<auxiliaryFunctionsMatrix<<std::endl;
//      std::cout<<"auxFunct row 4 = "<<auxiliaryFunctionsMatrix.row(4)<<std::endl;
//      std::cout<<"auxFunct row 5 = "<<auxiliaryFunctionsMatrix.row(5)<<std::endl;
//      std::cout<<"auxFunct row 6 = "<<auxiliaryFunctionsMatrix.row(6)<<std::endl;
//      std::cout<<"auxFunct row 7 = "<<auxiliaryFunctionsMatrix.row(7)<<std::endl;
//      std::cout<<"auxFunct row 8 = "<<auxiliaryFunctionsMatrix.row(8)<<std::endl;
//      std::cout<<"auxFunct row 9 = "<<auxiliaryFunctionsMatrix.row(9)<<std::endl;
//      std::cout<<"auxFunct row 27 = "<<auxiliaryFunctionsMatrix.row(27)<<std::endl;

    return auxiliaryFunctionsMatrix;

}






private:

    // The diferent celestial body constant parameters and polynomial coefficient parameter matrices

 double adiabeticIndex;                                   // gamma_a      adiabetic index
 double specificGasConstant;                        // Rstar    [m^2/(s^2*K)]    specific gas constant
 double standardGravitationalParameter;  // mu_M     [m^3/s^2]    standard gravitational parameter
 double rotationalVelocity;                         // rotational velocity of Mars  [rad/s]
 double primeMeridianAngle;                          // OmegaP   [rad]   relative angle between the prime meridian and the x-axis
 double inertialFrameTime;                            // t0       [s]    time between the start time and the time that the inertial frame was set
 double bodyReferenceRadius;                         // Rm    [m]       MOLA radius of Mars

 Eigen::MatrixXd temperaturePolyCoefficients; // PTn    temperature polynomial coefficients
 Eigen::MatrixXd temperatureAltitudeRanges;    // altitude range per section for the temperature-altitude curve [km MOLA]
 Eigen::VectorXd densityPolyCoefficients;         // Prho n density polynomial coefficients

     // The different vehicle constant parameters and polynomial coefficients

 double Thrust;                                                         // T   [N]  engine nominal thrust
 Eigen::MatrixXd thrustAzimuthMatrix;                           // psi_T [rad] thrust azimuth angles
 Eigen::MatrixXd thrustElevationMatrix;                         // epsilon_T [rad] thrust elevation angles
 double specificImpulse;                                        // Isp [s]    engine nominal specific impulse
 double referenceArea;                                           // S [m^2]  vehicle reference area
 Eigen::MatrixXd dragCoefficientPolyCoefficients;           // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
 Eigen::MatrixXd dragCoefficientMachRanges;                       // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient

    // Set functions

  double FlightPathAngle;         // Flight path angle in rad
 double HeadingAngle;            // Heading angle in rad



    // Additional in-class used variables


    int sectionT;                                   // This variable holds the "(section number -1)" for the temperature curve fit
    int powerT;                                     // This variable holds the section corresponding order for the temperature curve fit

    int sectionCD;                                  // This variable holds the "(section number -1)" for the drag coefficient curve fit

    Eigen::VectorXd auxiliaryEquationsVector;               // The vector containing the auxiliary equations xn
    Eigen::VectorXd auxiliaryDerivativesVector;             // The vector containing the auxiliary derivatives un
    Eigen::MatrixXd auxiliaryFunctionsMatrix;               // The matrix containing the auxiliary functions wn,m

};

#endif // AUXILIARY_H
