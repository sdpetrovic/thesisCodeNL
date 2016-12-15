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
 *      160719    S.D. Petrovic     Wrote everything in spherical coordinates
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

class Auxiliary

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
     * const Eigen::MatrixXd thrustAzimuthPolyCoefficients_;                 // P_psiT   these are the polynomial coefficients for the fit for the psiT curve
     * const Eigen::MatrixXd thrustElevationPolyCoefficients_;               // P_epsilonT   these are the polynomial coefficients for the fit for the epsilonT curve
     *
     *
     *
     */


    Auxiliary(const double adiabeticIndex_, const double specificGasConstant_, const double standardGravitationalParameter_, const double rotationalVelocity_, const double primeMeridianAngle_,
              const double inertialFrameTime_, const double bodyReferenceRadius_, const Eigen::MatrixXd temperaturePolyCoefficients_, const Eigen::MatrixXd temperatureAltitudeRanges_,
              const Eigen::VectorXd densityPolyCoefficients_, const double Thrust_, const Eigen::MatrixXd thrustAzimuthMatrix_, const Eigen::MatrixXd thrustElevationMatrix_, const double specificImpulse_,
              const double referenceArea_, const Eigen::MatrixXd dragCoefficientPolyCoefficients_, const Eigen::MatrixXd dragCoefficientMachRanges_, const Eigen::MatrixXd thrustAzimuthPolyCoefficients_,const Eigen::MatrixXd thrustElevationPolyCoefficients_){

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
         thrustAzimuthPolyCoefficients = thrustAzimuthPolyCoefficients_;            // P_psiT   these are the polynomial coefficients for the fit for the psiT curve
         thrustElevationPolyCoefficients = thrustElevationPolyCoefficients_;        // P_epsilonT   these are the polynomial coefficients for the fit for the epsilonT curve



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

    Eigen::VectorXd getAuxiliaryEquations( const tudat::basic_mathematics::Vector7d& aState, const double time, const Eigen::Vector3d& thrustAccelerationsBframe){
//        std::cout<<"verticalInertialFlightPathAngleSet eq 1 = "<<verticalInertialFlightPathAngleSet<<std::endl;
//        std::cout<<"verticalRotationalFlightPathAngleSet eq 1 = "<<verticalRotationalFlightPathAngleSet<<std::endl;

        auxiliaryEquationsVector = Eigen::VectorXd::Zero(35);       // Setting the complete vector and filling it with zeros for now



        // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry

        auxiliaryEquationsVector(11) = aState(2);              // x11 Longitude (tau)  [rad]
        auxiliaryEquationsVector(12) = aState(1);              // x12 Latitude (delta)  [rad]
        auxiliaryEquationsVector(13) = aState(5);              // x13 Azimuth angle (chi_G)    [rad
        auxiliaryEquationsVector(14) = aState(4);              // x14 Flight-path angle (gamma_G)  [rad]
        auxiliaryEquationsVector(15) = aState(3);              // x15 Ground velocity (V_G) (initially 0.0 [km/s])
        auxiliaryEquationsVector(16) = aState(0);              // x16 Radius (r)  [km]
        auxiliaryEquationsVector(7) = aState(6);              // x7

       //std::cout<<"Surely this works 1..."<<std::endl;

//        auxiliaryEquationsVector(8) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2)+
//                auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3) ;              // x8

//        auxiliaryEquationsVector(9) = pow(auxiliaryEquationsVector(8), 1.5);              // x9

//        auxiliaryEquationsVector(10) = rotationalVelocity*(inertialFrameTime+time)-primeMeridianAngle ;              // x10

//        const double inertialVelocity = sqrt(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(5)*auxiliaryEquationsVector(5)+auxiliaryEquationsVector(6)*auxiliaryEquationsVector(6)); //

//        // Account for rounding errors
//        if (abs(inertialVelocity*inertialVelocity+rotationalVelocity*rotationalVelocity*
//                (auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))+
//                2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1))) <= 1e-16){
//            auxiliaryEquationsVector(15) = 0.0;
//        }
//        else{
//        auxiliaryEquationsVector(15) = sqrt(inertialVelocity*inertialVelocity+rotationalVelocity*rotationalVelocity*
//                                            (auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))+
//                                            2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1)));               // x15
//}
//        /// Debug ///
//        std::cout<<"VI^2+Omega_M^2(x1^2+x2^2)+2*Omega_M*(x4*x2-x5*x1) = "<<inertialVelocity*inertialVelocity+rotationalVelocity*rotationalVelocity*
//                   (auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))+
//                   2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1))<<std::endl;

//        std::cout<<"x4*x2-x5*x1 = "<<auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1)<<std::endl;
//        std::cout<<"x1^2+x2^2 = "<<auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2)<<std::endl;
//        std::cout<<"VI^2+Omega_M^2(x1^2+x2^2) = "<<inertialVelocity*inertialVelocity+rotationalVelocity*rotationalVelocity*
//                   (auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2))<<std::endl;
//        std::cout<<"2*Omega_M*(x4*x2-x5*x1) = "<<2.0*rotationalVelocity*(auxiliaryEquationsVector(4)*auxiliaryEquationsVector(2)-auxiliaryEquationsVector(5)*auxiliaryEquationsVector(1))<<std::endl;
//        /// Debug ///


//        auxiliaryEquationsVector(16) = sqrt(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(1)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(2)+auxiliaryEquationsVector(3)*auxiliaryEquationsVector(3));               // x16

        /// Debug ///
//        std::cout<<"x16 = "<<auxiliaryEquationsVector(16)<<std::endl;
        /// Debug ///

        // Avoid cosine rounding errors
//        std::cout<<"cos(pi/2) = "<<cos(tudat::mathematical_constants::LONG_PI/2.0)<<std::endl;
//        double cx10;

//        if (abs(cos(auxiliaryEquationsVector(10)))<6.2e-17){
//            cx10 = 0;
//        }
//        else {
//            cx10 = cos(auxiliaryEquationsVector(10));
//        }

        // Same for sin
//        std::cout<<"sin(pi) = "<<sin(tudat::mathematical_constants::LONG_PI)<<std::endl;

//       double sx10;

//        if (abs(sin(auxiliaryEquationsVector(10)))<1.2e-16){

//            sx10 = 0;

//        }
//        else {
//            sx10 = sin(auxiliaryEquationsVector(10));
//        }

//        const double rotationalXposition = cx10*auxiliaryEquationsVector(1)+sx10*auxiliaryEquationsVector(2); // x_R
//        const double rotationalYposition = -sx10*auxiliaryEquationsVector(1)+cx10*auxiliaryEquationsVector(2); // y_R

//        auxiliaryEquationsVector(11) = atan2(rotationalYposition,rotationalXposition);               // x11

//         auxiliaryEquationsVector(12) = asin(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(16));               // x12

         // Avoid cosine rounding errors
//                  double cx11;

////                  std::cout<<"x11 = "<<auxiliaryEquationsVector(11)<<std::endl;

//         if (abs(cos(auxiliaryEquationsVector(11)))<6.2e-17){
//             cx11 = 0;
//         }
//         else {
//             cx11 = cos(auxiliaryEquationsVector(11));

////             std::cout<<"x11 has been set"<<std::endl;
//         }

//         double cx12;

//         if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17){
//             cx12 = 0;
//         }
//         else {
//             cx12 = cos(auxiliaryEquationsVector(12));
//         }



         // Same for sin

//         double sx11;

//          if (abs(sin(auxiliaryEquationsVector(11)))<1.2e-16){

//              sx11 = 0;

//          }
//          else {
//              sx11 = sin(auxiliaryEquationsVector(11));
//          }

//         double sx12;

//         if (abs(sin(auxiliaryEquationsVector(12)))<6.2e-17){
//             sx12 = 0;
//         }
//         else {
//             sx12 = sin(auxiliaryEquationsVector(12));
//         }


//        const double verticalXvelocity = (auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2))*(sx12*sx10*sx11-
//                                                                                                                       cx10*cx11*sx12)+
//                (auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1))*(-cx11*sx12*sx10-
//                                                                                              cx10*sx12*sx11)+
//                auxiliaryEquationsVector(6)*cx12; // Vx_V

//        const double verticalYvelocity = (auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2))*(-cx11*sx10-
//                                                                                                                       cx10*sx11)+
//                (auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1))*(cx10*cx11-
//                                                                                              sx10*sx11); // Vy_V

//        const double verticalZvelocity = (auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2))*(cx12*sx10*sx11-
//                                                                                                                       cx10*cx11*cx12)+
//                (auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1))*(-cx11*cx12*sx10-
//                                                                                              cx10*cx12*sx11)-
//                auxiliaryEquationsVector(6)*sx12; // Vz_V

        /// Debug ///
//        std::cout<<"x10 = "<<auxiliaryEquationsVector(10)<<std::endl;
//        std::cout<<"x11 = "<<auxiliaryEquationsVector(11)<<std::endl;
//        std::cout<<"cx10 = "<<cx10<<std::endl;
//        std::cout<<"sx10 = "<<sx10<<std::endl;
//        std::cout<<"cx11 = "<<cx11<<std::endl;
//        std::cout<<"sx11 = "<<sx11<<std::endl;
//        std::cout<<"x4+Omega_M*x2 = "<<(auxiliaryEquationsVector(4)+rotationalVelocity*auxiliaryEquationsVector(2))<<std::endl;
//        std::cout<<"x5 - Omega_M*x1 = "<<(auxiliaryEquationsVector(5)-rotationalVelocity*auxiliaryEquationsVector(1))<<std::endl;
//        std::cout<<"y_R = "<<rotationalYposition<<std::endl;
//        std::cout<<"x_R = "<<rotationalXposition<<std::endl;
//        std::cout<<"verticalXvelocity = "<<verticalXvelocity<<std::endl;
//        std::cout<<"verticalYvelocity = "<<verticalYvelocity<<std::endl;
//        std::cout<<"verticalZvelocity = "<<verticalZvelocity<<std::endl;
//        std::cout<<"V_G = "<<auxiliaryEquationsVector(15)<<std::endl;
//        std::cout<<"V_G - verticalZvelocity = "<<auxiliaryEquationsVector(15)+verticalZvelocity<<std::endl;
//        std::cout<<"verticalZvelocity/V_G = "<<verticalZvelocity/auxiliaryEquationsVector(15)<<std::endl;
//        std::cout<<"verticalZvelocity/V_G+1 = "<<verticalZvelocity/auxiliaryEquationsVector(15)+1<<std::endl;
//        std::cout<<"FlightPathAngle = "<<FlightPathAngle<<std::endl;
//        std::cout<<"HeadingAngle = "<<HeadingAngle<<std::endl;
        /// Debug ///

//        // Check if the angles have been set and then use them as initial condition
//        if (FlightPathAngle != 1000 && HeadingAngle != 1000 && time == 0.0){
//            auxiliaryEquationsVector(13) = HeadingAngle;
//            auxiliaryEquationsVector(14) = FlightPathAngle;
//        }
//        else {

//        auxiliaryEquationsVector(13) = atan2(verticalYvelocity,verticalXvelocity);               // x13

////        /// Debug ///
//        std::cout<<"x13 first = "<<auxiliaryEquationsVector(13)<<std::endl;
////        std::cout<<"atan2(-1,0) = "<<atan2(-1,0)<<std::endl;
////        /// Debug ///

//        // Take care of the 0 m/s velocity singularity
//        if (auxiliaryEquationsVector(15) == 0){
//            auxiliaryEquationsVector(14) = tudat::mathematical_constants::LONG_PI/2.0;
//        }
//        else if (verticalZvelocity/auxiliaryEquationsVector(15) >= 1.0 || verticalZvelocity/auxiliaryEquationsVector(15)-1 >= -1E-15){        // Compensate for rounding errors
//            auxiliaryEquationsVector(14) = -asin(1.0);
//        }
//        else if (verticalZvelocity/auxiliaryEquationsVector(15) <= -1.0 || verticalZvelocity/auxiliaryEquationsVector(15)+1 <= 1E-15){       // Compensate for rounding errors
//            auxiliaryEquationsVector(14) = -asin(-1.0);
//        }
//        else {
//        auxiliaryEquationsVector(14) = -asin(verticalZvelocity/auxiliaryEquationsVector(15));               // x14
//} // End of check for set angles
        /// Debug ///

//        std::cout<<"x14 = "<<auxiliaryEquationsVector(14)<<std::endl;

        /// Debug ///


//        std::cout<<"verticalZvelocity/auxiliaryEquationsVector(15) -1 = "<<verticalZvelocity/auxiliaryEquationsVector(15) -1 <<std::endl;
//}





        // Please note that the altitude h (x31) is expressed in km MOLA (which is also the input for the density and temperature curves!)
//        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(20)-bodyReferenceRadius)/1000;              // x31 [km]!!!
//        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(20)*1e6-bodyReferenceRadius*1e6)/1e6;              // x31 [km]!!!
        auxiliaryEquationsVector(31) = (auxiliaryEquationsVector(16)-bodyReferenceRadius);              // x31 [km]!!!

        // Determining the Thrust azimuth (psiT) and the Thrust elevation (epsilonT) angles based on the altitude
                // Determine the proper azimuth value for the current altitude section
                int sectionThrustAz = 0;    // Set the current azimuth value to the default first section
                for (int i = 0; i < thrustAzimuthMatrix.rows();i++){
                    if (thrustAzimuthMatrix(i,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < thrustAzimuthMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                        sectionThrustAz = i;
//                        std::cout<<"sectionThrustAz = "<<sectionThrustAz+1<<std::endl;
                    }
                }

                const double thrustAzimuth = thrustAzimuthMatrix(sectionThrustAz,2); // Set the thrust azimuth to the current azimuth corresponding to the current altitude section

                // Determine the proper elevation value for the current altitude section
                int sectionThrustEl = 0;    // Set the current elevation value to the default first section
                for (int i = 0; i < thrustElevationMatrix.rows();i++){
                    if (thrustElevationMatrix(i,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < thrustElevationMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                        sectionThrustEl = i;
                    }
                }

                const double thrustElevation = thrustElevationMatrix(sectionThrustEl,2); // Set the thrust elevation to the current elevation corresponding to the current altitude section

//                const double thrustAzimuth = 0.0; // Test
//                const double thrustElevation = 0.0; // Test


        // Computing the polynomial fit using the altitude and fit parameters for density
        for (int i = 0; i < 10+1;i++) {

        auxiliaryEquationsVector(30) += pow(auxiliaryEquationsVector(31),i)*densityPolyCoefficients(i);              // x30
};

        /// Debug ///

//        double testLNrho = 0.0;

//        for (int i = 0; i < 10+1;i++) {

//        testLNrho += pow(380.0,i)*densityPolyCoefficients(i);              // test
//};

//        std::cout<<"auxiliaryEquationsVector(30) (LNrho) = "<<auxiliaryEquationsVector(30)<<std::endl;
//        std::cout<<"testLNrho = "<<testLNrho<<std::endl;

        /// Debug ///

        // Determine which section of the temperature curve needs to be used and what the corresponding order is
        // Also, because a computer is less than perfect, a small correction is made to the lower bound of the first section to make sure that the initial altitude is still valid

        if ((temperatureAltitudeRanges(0,0)-0.000000000001) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(0,1)){

        sectionT = 0;
        powerT = 1;

        }
        else if (temperatureAltitudeRanges(1,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(1,1)){

        sectionT = 1;
        powerT = 3;

        }
        else if (temperatureAltitudeRanges(2,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(2,1)){

        sectionT = 2;
        powerT = 6;

        }
        else if (temperatureAltitudeRanges(3,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(3,1)){

            sectionT = 3;
            powerT = 8;
        }
        else if (temperatureAltitudeRanges(4,0) <= auxiliaryEquationsVector(31)){

            sectionT = 4;
            powerT = 0;
        }
        else {


            std::cerr<<"The current altitude: "<<auxiliaryEquationsVector(31)<<" [km MOLA] is not a valid altitude (lower than the lowest reference altitude)"<<std::endl;

                       sectionT = 0;
                        powerT = 1;

        };
        //std::cout<<"Surely this works 4..."<<std::endl;

        // Computing the polynomial fit using the altitude and fit parameters for temperature
        for (int i=0; i < powerT+1;i++){

        auxiliaryEquationsVector(34) += pow(auxiliaryEquationsVector(31),i)*temperaturePolyCoefficients(sectionT,i);              // x34

//        std::cout<<"x34 interval "<<i<<" = "<<auxiliaryEquationsVector(34)<<std::endl;

};


        auxiliaryEquationsVector(28) = exp(auxiliaryEquationsVector(30));              // x28

        if (referenceArea == 0.0){ // In case the density equation does weird things
            auxiliaryEquationsVector(28) = 0.0;
        }

        auxiliaryEquationsVector(33) = sqrt(adiabeticIndex*specificGasConstant*auxiliaryEquationsVector(34));              // x33

//std::cout<<"Surely this works 5..."<<std::endl;

        auxiliaryEquationsVector(32) = auxiliaryEquationsVector(15)/auxiliaryEquationsVector(33);              // x32


        // Determine which section of the drag coefficient curve needs to be used

        /// Debug ///
//        std::cout<<"dragCoefficientMachRanges = "<<dragCoefficientMachRanges<<std::endl;
        /// Debug ///

        for (int i=0; i < 5+1; i++){

            if (dragCoefficientMachRanges(i,0) <= auxiliaryEquationsVector(32) && auxiliaryEquationsVector(32) < dragCoefficientMachRanges(i,1)){

                sectionCD = i;

//                std::cout<<"sectionCD is set equal to "<<sectionCD<<std::endl;


            }


        };

//std::cout<<"Surely this works 6..."<<std::endl;

/// Debug ///
//std::cout<<"sectionCD = "<<sectionCD<<std::endl;
//std::cout<<"x32 = "<<auxiliaryEquationsVector(32)<<std::endl;
/// Debug ///

        auxiliaryEquationsVector(29) = dragCoefficientPolyCoefficients(sectionCD,1)*auxiliaryEquationsVector(32)+dragCoefficientPolyCoefficients(sectionCD,0);              // x29

//std::cout<<"Surely this works 7..."<<std::endl;

        auxiliaryEquationsVector(27) = 0.5*referenceArea*auxiliaryEquationsVector(28)*auxiliaryEquationsVector(15)*auxiliaryEquationsVector(15)*auxiliaryEquationsVector(29);              // x27

        if (referenceArea == 0.0){ // In case the density equation does weird things
            auxiliaryEquationsVector(27) = 0.0;
        }

//        // Accounting for altitude higher than 320 km MOLA: assumption, density can be neglected
//        if (auxiliaryEquationsVector(31)>=320.0){
//            auxiliaryEquationsVector(27) = 0.0;
//            std::cout<<"It's done, it's gone!"<<std::endl;
//        }

//std::cout<<"Surely this works 8..."<<std::endl;

//        auxiliaryEquationsVector(0) = thrustAccelerationsBframe(0)-(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7));              // w4,2
        auxiliaryEquationsVector(0) = (Thrust*cos(thrustAzimuth)*cos(thrustElevation)/auxiliaryEquationsVector(7))-(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7));              // w4,2

//        /// Debug ///
//        std::cout<<"auxiliaryEquationsVector(0) = "<<auxiliaryEquationsVector(0)<<std::endl;
//        std::cout<<"Thrust*cos(thrustAzimuth)*cos(thrustElevation)/auxiliaryEquationsVector(7) = "<<Thrust*cos(thrustAzimuth)*cos(thrustElevation)/auxiliaryEquationsVector(7)<<std::endl;
//        std::cout<<"-(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7)) = "<<-(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7))<<std::endl;
//        std::cout<<"auxiliaryEquationsVector(27) (D) = "<<auxiliaryEquationsVector(27)<<std::endl;
//        std::cout<<"referenceArea = "<<referenceArea<<std::endl;
//        std::cout<<"auxiliaryEquationsVector(28) (rho) = "<<auxiliaryEquationsVector(28)<<std::endl;
//        std::cout<<"auxiliaryEquationsVector(15) (V_G) = "<<auxiliaryEquationsVector(15)<<std::endl;
//        std::cout<<"auxiliaryEquationsVector(29) (C_D) = "<<auxiliaryEquationsVector(29)<<std::endl;
//        std::cout<<"auxiliaryEquationsVector(31) (h) = "<<auxiliaryEquationsVector(31)<<std::endl;
//        std::cout<<"auxiliaryEquationsVector(7) = "<<auxiliaryEquationsVector(7)<<std::endl;

//        std::cout<<"thrustAccelerationsBframe(0) = "<<thrustAccelerationsBframe(0)<<std::endl;
//        std::cout<<"(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7)) = "<<(auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7))<<std::endl;
//        /// Debug ///



// Set vertical ascent to false again
//        verticalInertialFlightPathAngleSet = false;
//        verticalInertialFlightPathAngleSet = NULL;

//        verticalRotationalFlightPathAngleSet = false;
//        verticalRotationalFlightPathAngleSet = NULL;



// auxiliaryEquationsVector() = ;               // x


//std::cout<<"Surely this works end..."<<std::endl;


       return auxiliaryEquationsVector;
    }

//////////////////////////////////////////////// Auxiliary Derivatives //////////////////////////////////////////////////////////////////////

    Eigen::VectorXd getAuxiliaryDerivatives( const tudat::basic_mathematics::Vector7d& aState, const double time, const Eigen::Vector3d& thrustAccelerationsBframe, const Eigen::VectorXd& auxiliaryEquationsVector){

//        std::cout<<"verticalInertialFlightPathAngleSet der 1 = "<<verticalInertialFlightPathAngleSet<<std::endl;
//        std::cout<<"verticalRotationalFlightPathAngleSet der 1 = "<<verticalRotationalFlightPathAngleSet<<std::endl;

    auxiliaryDerivativesVector = Eigen::VectorXd::Zero(35);       // Setting the complete vector and filling it with zeros for now

    // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry
    // Which in this case means that the first entry of the vector is 0 and is not used.

//    auxiliaryDerivativesVector(10) = 0;
//    auxiliaryDerivativesVector(10) = rotationalVelocity;                // u10


//    auxiliaryDerivativesVector(1) = auxiliaryEquationsVector(4);                // u1

//    auxiliaryDerivativesVector(2) = auxiliaryEquationsVector(5);                // u2

//    auxiliaryDerivativesVector(3) = auxiliaryEquationsVector(6);                // u3



    // Avoid cosine rounding errors
//    double cx10x11;

//    if (abs(cos(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11)))<6.2e-17){
//        cx10x11 = 0;
//    }
//    else {
//        cx10x11 = cos(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11));
//    }


    // Determining the Thrust azimuth (psiT) and the Thrust elevation (epsilonT) angles based on the altitude
            // Determine the proper azimuth value for the current altitude section
            int sectionThrustAz = 0;    // Set the current azimuth value to the default first section
            for (int i = 0; i < thrustAzimuthMatrix.rows();i++){
                if (thrustAzimuthMatrix(i,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < thrustAzimuthMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                    sectionThrustAz = i;
//                    std::cout<<"sectionThrustAz = "<<sectionThrustAz+1<<std::endl;
                }
            }

                const double thrustAzimuth = thrustAzimuthMatrix(sectionThrustAz,2); // Set the thrust azimuth to the current azimuth corresponding to the current altitude section

            // Determine the proper elevation value for the current altitude section
            int sectionThrustEl = 0;    // Set the current elevation value to the default first section
            for (int i = 0; i < thrustElevationMatrix.rows();i++){
                if (thrustElevationMatrix(i,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < thrustElevationMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                    sectionThrustEl = i;
                }
            }

                const double thrustElevation = thrustElevationMatrix(sectionThrustEl,2); // Set the thrust elevation to the current elevation corresponding to the current altitude section

//            const double thrustAzimuth = 0.0; // Test
//            const double thrustElevation = 0.0; // Test


    double cx12;

    if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17){
        cx12 = 0;
    }
    else {
        cx12 = cos(auxiliaryEquationsVector(12));
    }

    double cx13;

    if (abs(cos(auxiliaryEquationsVector(13)))<6.2e-17){
        cx13 = 0;
    }
    else {
        cx13 = cos(auxiliaryEquationsVector(13));
    }

    double cx14;

    if (abs(cos(auxiliaryEquationsVector(14)))<6.2e-17){
        cx14 = 0;
//        std::cout<<"cx14 has been set to 0"<<std::endl;
    }
    else {
        cx14 = cos(auxiliaryEquationsVector(14));
    }

    // Same for sin

//   double sx10x11;

//    if (abs(sin(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11)))<1.2e-16){ // Using x10 and x11

//        sx10x11 = 0;

//    }
//    else {
//        sx10x11 = sin(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11));
//    }

    double sx12;

    if (abs(sin(auxiliaryEquationsVector(12)))<6.2e-17){
        sx12 = 0;
    }
    else {
        sx12 = sin(auxiliaryEquationsVector(12));
    }

    double sx13;

    if (abs(sin(auxiliaryEquationsVector(13)))<6.2e-17){
        sx13 = 0;
    }
    else {
        sx13 = sin(auxiliaryEquationsVector(13));
    }

    double sx14;

    if (abs(sin(auxiliaryEquationsVector(14)))<6.2e-17){
        sx14 = 0;
//        std::cout<<"cx14 has been set to 0"<<std::endl;
    }
    else {
        sx14 = sin(auxiliaryEquationsVector(14));
    }


//    auxiliaryDerivativesVector(4) = -standardGravitationalParameter*(auxiliaryEquationsVector(1)/auxiliaryEquationsVector(9))+auxiliaryEquationsVector(0)*
//            (cx10x11*(-sx12*cx13*cx14+cx12*sx14)-
//                 sx10x11*sx13*cx14)+
//            thrustAccelerationsBframe(1)*(cx10x11*sx12*sx13-sx10x11*cx13)+
//            thrustAccelerationsBframe(2)*(cx10x11*(-sx12*cx13*sx14-
//                                                                                                          cx12*cx14)-
//                                          sx10x11*sx13*sx14);                // u4




//    auxiliaryDerivativesVector(5) = -standardGravitationalParameter*(auxiliaryEquationsVector(2)/auxiliaryEquationsVector(9))+auxiliaryEquationsVector(0)*
//            (sx10x11*(-sx12*cx13*cx14+cx12*sx14)+
//                 cx10x11*sx13*cx14)+
//            thrustAccelerationsBframe(1)*(sx10x11*sx12*sx13+cx10x11*cx13)+
//            thrustAccelerationsBframe(2)*(sx10x11*(-sx12*cx13*sx14-
//                                                                                                          cx12*cx14)+
//                                          cx10x11*sx13*sx14);                // u5

//    auxiliaryDerivativesVector(6) = -standardGravitationalParameter*(auxiliaryEquationsVector(3)/auxiliaryEquationsVector(9))+auxiliaryEquationsVector(0)*
//            (cx12*cx13*cx14+sx12*sx14)-
//            thrustAccelerationsBframe(1)*cx12*sx13+
//            thrustAccelerationsBframe(2)*(cx12*cx13*sx14-sx12*cx14);                // u6



    auxiliaryDerivativesVector(7) = -(Thrust/(tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION*specificImpulse));                // u7

//    auxiliaryDerivativesVector(8) = 2.0*(auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4)+auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5)+auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6));                // u8


    // Account for singularities at the poles
    if (auxiliaryEquationsVector(12) == tudat::mathematical_constants::LONG_PI/2.0 || auxiliaryEquationsVector(12) == -tudat::mathematical_constants::LONG_PI/2.0 || cx12 == 0.0){
        auxiliaryDerivativesVector(11) = 0;
    }
    else {
    auxiliaryDerivativesVector(11) = (auxiliaryEquationsVector(15)*sx13*cx14)/(auxiliaryEquationsVector(16)*cx12);                // u11
}


    auxiliaryDerivativesVector(12) = (auxiliaryEquationsVector(15)*cx13*cx14)/auxiliaryEquationsVector(16);                // u12


//    /// Debug ///
//    std::cout<<"x13 = "<<auxiliaryEquationsVector(13)<<std::endl;
//    std::cout<<"x14 = "<<auxiliaryEquationsVector(14)<<std::endl;
//    std::cout<<"x14 - pi/2 = "<<auxiliaryEquationsVector(14)-tudat::mathematical_constants::LONG_PI/2.0<<std::endl;
//    std::cout<<"x12 = "<<auxiliaryEquationsVector(12)<<std::endl;
//    std::cout<<"x12 - pi/2 = "<<auxiliaryEquationsVector(12)-tudat::mathematical_constants::LONG_PI/2.0<<std::endl;
//    /// Debug ///


    // Account for singularities at vertical ascent and zero velocity
    if (auxiliaryEquationsVector(14) == tudat::mathematical_constants::LONG_PI/2.0 || auxiliaryEquationsVector(14) == -tudat::mathematical_constants::LONG_PI/2.0 || cx14 == 0.0){ // Flight-path angle
        if (abs(auxiliaryEquationsVector(12)-tudat::mathematical_constants::LONG_PI/2.0) <= 1e-16 || abs(auxiliaryEquationsVector(12)+tudat::mathematical_constants::LONG_PI/2.0) <= 1e-16 ){
            auxiliaryDerivativesVector(13) = 0.0;
        }
        else {
        auxiliaryDerivativesVector(13) = (auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cx14*tan(auxiliaryEquationsVector(12))*sx13;
        }

    }
    else if (auxiliaryEquationsVector(15) == 0){ // Velocity
        auxiliaryDerivativesVector(13) = 2.0*(rotationalVelocity/cx14)*(sx12*cx14-
                                                                                                     cx12*sx14*cx13);
    }
    else if (abs(auxiliaryEquationsVector(12)-tudat::mathematical_constants::LONG_PI/2.0) <= 1e-16 || abs(auxiliaryEquationsVector(12)+tudat::mathematical_constants::LONG_PI/2.0) <= 1e-16 ){
       auxiliaryDerivativesVector(13) = 2.0*(rotationalVelocity/cx14)*(sx12*cx14-
                                                                                                         cx12*sx14*cx13)+
                    (rotationalVelocity*rotationalVelocity/(auxiliaryEquationsVector(15)*cx14))*auxiliaryEquationsVector(16)*cx12*sx12*sx13+
//                    (thrustAccelerationsBframe(1)/(auxiliaryEquationsVector(15)*cx14));
               ((Thrust*sin(thrustAzimuth)*cos(thrustElevation))/(auxiliaryEquationsVector(15)*auxiliaryEquationsVector(7)*cx14));
    }
    else {
    auxiliaryDerivativesVector(13) = 2.0*(rotationalVelocity/cx14)*(sx12*cx14-
                                                                                                 cx12*sx14*cx13)+
            (auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cx14*tan(auxiliaryEquationsVector(12))*sx13+
            (rotationalVelocity*rotationalVelocity/(auxiliaryEquationsVector(15)*cx14))*auxiliaryEquationsVector(16)*cx12*sx12*sx13+
//            (thrustAccelerationsBframe(1)/(auxiliaryEquationsVector(15)*cx14));                // u13
            ((Thrust*sin(thrustAzimuth)*cos(thrustElevation))/(auxiliaryEquationsVector(15)*auxiliaryEquationsVector(7)*cx14));               // u13
}

    // Account for singularities at zero velocity
    if (auxiliaryEquationsVector(15) == 0.0){
        auxiliaryDerivativesVector(14) = 2.0*rotationalVelocity*cx12*sx13;
//        std::cout<<"u14 goes here 1"<<std::endl;
    }
    else {
//        std::cout<<"u14 goes here 2"<<std::endl;
//    auxiliaryDerivativesVector(14) = 2.0*rotationalVelocity*cx12*sx13+(auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cx14+
//            (rotationalVelocity*rotationalVelocity/auxiliaryEquationsVector(15))*auxiliaryEquationsVector(16)*cx12*(cx12*cx14+sx14*sx12*cx13)
//            -thrustAccelerationsBframe(2)/auxiliaryEquationsVector(15)-standardGravitationalParameter*cx14/(auxiliaryEquationsVector(15)*auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16));                // u14

    auxiliaryDerivativesVector(14) = 2.0*rotationalVelocity*cx12*sx13+(auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cx14+
            ((rotationalVelocity*rotationalVelocity)*auxiliaryEquationsVector(16)*cx12*(cx12*cx14+sx14*sx12*cx13)+
//            -thrustAccelerationsBframe(2)-standardGravitationalParameter*cx14/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16)))/auxiliaryEquationsVector(15);                // u14
             (Thrust*sin(thrustElevation))/(auxiliaryEquationsVector(7))-standardGravitationalParameter*cx14/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16)))/auxiliaryEquationsVector(15);                // u14

//    auxiliaryDerivativesVector(14) = 2.0*rotationalVelocity*cx12*sx13+(auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cx14+
//            ((-standardGravitationalParameter*cx14)/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))-thrustAccelerationsBframe(2)+
//             (rotationalVelocity*rotationalVelocity)*auxiliaryEquationsVector(16)*cx12*(cx12*cx14+sx14*sx12*cx13))/(auxiliaryEquationsVector(15)); // u14

    }
    /// Debug ///
//    std::cout<<"x10 = "<<auxiliaryEquationsVector(10)<<std::endl;
//    std::cout<<"x11 = "<<auxiliaryEquationsVector(11)<<std::endl;
//    std::cout<<"tau (using lambda) = "<<atan2(auxiliaryEquationsVector(2),auxiliaryEquationsVector(1))-auxiliaryEquationsVector(10)<<std::endl;
//    std::cout<<"x12 = "<<auxiliaryEquationsVector(12)<<std::endl;
//    std::cout<<"u14 = "<<auxiliaryDerivativesVector(14)<<std::endl;
//    std::cout<<"u14 = "<<2.0*rotationalVelocity*cx12*sx13+(auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cx14+
//               (rotationalVelocity*rotationalVelocity/auxiliaryEquationsVector(15))*auxiliaryEquationsVector(16)*cx12*(cx12*cx14+
//                                                                                                                                                    sx14*sx12*cx13)+
//               -thrustAccelerationsBframe(2)/auxiliaryEquationsVector(15)-standardGravitationalParameter*cx14/(auxiliaryEquationsVector(15)*auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))<<std::endl;
//    std::cout<<"u14 part 1 = "<<2.0*rotationalVelocity*cx12*sx13<<std::endl;
//    std::cout<<"u14 part 2 = "<<(auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cx14<<std::endl;
//    std::cout<<"u14 part 3 = "<<(rotationalVelocity*rotationalVelocity/auxiliaryEquationsVector(15))*auxiliaryEquationsVector(16)*cx12*(cx12*cx14+sx14*sx12*cx13)<<std::endl;
//    std::cout<<"u14 part 4 = "<<thrustAccelerationsBframe(2)/auxiliaryEquationsVector(15)<<std::endl;
//    std::cout<<"u14 part 5 = "<<standardGravitationalParameter*cx14/(auxiliaryEquationsVector(15)*auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))<<std::endl;
//    std::cout<<"w14,0 = "<<-thrustAccelerationsBframe(2)<<std::endl;
//    std::cout<<"w14,1 = "<<cx12*cx14+sx14*sx12*cx13<<std::endl;
//    std::cout<<"w14,2 = "<<-standardGravitationalParameter*cx14<<std::endl;
//    std::cout<<"w14,3 = "<<auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16)<<std::endl;
//    std::cout<<"w14,4 = "<<(-standardGravitationalParameter*cx14)/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))<<std::endl;
//    std::cout<<"w14,5 = "<<(-standardGravitationalParameter*cx14)/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))-thrustAccelerationsBframe(2)+(rotationalVelocity*rotationalVelocity)*auxiliaryEquationsVector(16)*cx12*(cx12*cx14+sx14*sx12*cx13)<<std::endl;
//    std::cout<<"w14,6 = "<<((-standardGravitationalParameter*cx14)/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))-thrustAccelerationsBframe(2)+(rotationalVelocity*rotationalVelocity)*auxiliaryEquationsVector(16)*cx12*(cx12*cx14+sx14*sx12*cx13))/(auxiliaryEquationsVector(15))<<std::endl;
//    std::cout<<"w14,7 = "<<2.0*rotationalVelocity*cx12*sx13+(auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cx14+((-standardGravitationalParameter*cx14)/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))-thrustAccelerationsBframe(2)+(rotationalVelocity*rotationalVelocity)*auxiliaryEquationsVector(16)*cx12*(cx12*cx14+sx14*sx12*cx13))/(auxiliaryEquationsVector(15))<<std::endl;


    /// Debug ///

    auxiliaryDerivativesVector(15) = rotationalVelocity*rotationalVelocity*auxiliaryEquationsVector(16)*cx12*
            (sx14*cx12-cx14*sx12*cx13)+
            auxiliaryEquationsVector(0)-standardGravitationalParameter*sx14/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16));                // u15

    /// Debug ///
//    std::cout<<"u15 = "<<auxiliaryDerivativesVector(15)<<std::endl;
//    std::cout<<"auxiliaryEquationsVector(0)-standardGravitationalParameter*sx14/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16)) = "<<auxiliaryEquationsVector(0)-standardGravitationalParameter*sx14/(auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))<<std::endl;
//    std::cout<<"auxiliaryEquationsVector(0) = "<<auxiliaryEquationsVector(0)<<std::endl;
    //    std::cout<<"u14 part 1 = "<<2.0*rotationalVelocity*cx12*sx13<<std::endl;
//    std::cout<<"u14 part 2 = "<<(auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cx14<<std::endl;
//    std::cout<<"u14 part 3 = "<<(rotationalVelocity*rotationalVelocity/auxiliaryEquationsVector(15))*auxiliaryEquationsVector(16)*cx12*(cx12*cx14+
//                                                                                                                                                                     sx14*sx12*cx13)<<std::endl;
//    std::cout<<"u14 part 4 = "<<-thrustAccelerationsBframe(2)/auxiliaryEquationsVector(15)<<std::endl;
//    std::cout<<"u14 part 5 = "<<standardGravitationalParameter*cx14/(auxiliaryEquationsVector(15)*auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))<<std::endl;

//    std::cout<<"u14 part2.1 = "<<(auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))<<std::endl;
//    std::cout<<"u14 part2.2 = "<<cx14<<std::endl;
               //    std::cout<<"u12 = "<<(auxiliaryEquationsVector(15)*cx13*cos(auxiliaryEquationsVector(14)))/auxiliaryEquationsVector(16)<<std::endl;
//    std::cout<<"x16 = "<<auxiliaryEquationsVector(16)<<std::endl;
//    std::cout<<"x15*cos(x13)*cos(x14) = "<<auxiliaryEquationsVector(15)*cx13*cos(auxiliaryEquationsVector(14))<<std::endl;
//    std::cout<<"x14 = "<<auxiliaryEquationsVector(14)<<std::endl;
//    std::cout<<"x14 - pi/2 = "<<auxiliaryEquationsVector(14)-tudat::mathematical_constants::LONG_PI/2.0<<std::endl;
//    std::cout<<"u14 = "<<2.0*rotationalVelocity*cx12*sx13+(auxiliaryEquationsVector(15)/auxiliaryEquationsVector(16))*cos(auxiliaryEquationsVector(14))+
//               (rotationalVelocity*rotationalVelocity/auxiliaryEquationsVector(15))*auxiliaryEquationsVector(16)*cx12*(cx12*cos(auxiliaryEquationsVector(14))+
//                                                                                                                                                    sin(auxiliaryEquationsVector(14))*sin(auxiliaryEquationsVector(12))*cx13)+
//               -thrustAccelerationsBframe(2)/auxiliaryEquationsVector(15)-standardGravitationalParameter*cos(auxiliaryEquationsVector(14))/(auxiliaryEquationsVector(15)*auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16))<<std::endl;

    /// Debug ///

    auxiliaryDerivativesVector(16) = auxiliaryEquationsVector(15)*sx14;                // u16


    /// Debug ///

//    std::cout<<"x14 = "<<auxiliaryEquationsVector(14)<<std::endl;
//    std::cout<<"x15 = "<<auxiliaryEquationsVector(15)<<std::endl;
//    std::cout<<"x16 = "<<auxiliaryEquationsVector(16)<<std::endl;
//    std::cout<<"u14 = "<<auxiliaryDerivativesVector(14)<<std::endl;
//    std::cout<<"u15 = "<<auxiliaryDerivativesVector(15)<<std::endl;
//    std::cout<<"u16 = "<<auxiliaryDerivativesVector(16)<<std::endl;

    /// Debug ///

//    auxiliaryDerivativesVector(9) = 1.5*((auxiliaryEquationsVector(9)*auxiliaryDerivativesVector(8))/auxiliaryEquationsVector(8));                // u9

    auxiliaryDerivativesVector(31) = auxiliaryDerivativesVector(16);                    // u31


    // Computing the polynomial fit derivative using the altitude and fit parameters for density
    for (int i = 1; i < 10+1;i++) {

    auxiliaryDerivativesVector(30) += auxiliaryDerivativesVector(31)*i*pow(auxiliaryEquationsVector(31),i-1)*densityPolyCoefficients(i);              // u30
};


    // Determine which section of the temperature curve needs to be used and what the corresponding order is
    // Also, because a computer is less than perfect, a small correction is made to the lower bound of the first section to make sure that the initial altitude is still valid

    if (temperatureAltitudeRanges(0,0)-0.000000000001 <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(0,1)){

    sectionT = 0;
    powerT = 1;

    }
    else if (temperatureAltitudeRanges(1,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(1,1)){

    sectionT = 1;
    powerT = 2;

    }
    else if (temperatureAltitudeRanges(2,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(2,1)){

    sectionT = 2;
    powerT = 6;

    }
    else if (temperatureAltitudeRanges(3,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < temperatureAltitudeRanges(3,1)){

        sectionT = 3;
        powerT = 8;
    }
    else if (temperatureAltitudeRanges(4,0) <= auxiliaryEquationsVector(31)){

        sectionT = 4;
        powerT = 0;
    }
    else {


        std::cerr<<"The current altitude: "<<auxiliaryEquationsVector(31)<<" [km MOLA] is not a valid altitude (lower than the lowest reference altitude)"<<std::endl;

                   sectionT = 0;
                    powerT = 1;

    };

    // Computing the polynomial fit using the altitude and fit parameters for temperature
    for (int i=1; i < powerT+1;i++){

    auxiliaryDerivativesVector(34) += auxiliaryDerivativesVector(31)*i*pow(auxiliaryEquationsVector(31),i-1)*temperaturePolyCoefficients(sectionT,i);              // u34

//    std::cout<<"i = "<<i<<std::endl;
//    std::cout<<"auxiliaryDerivativesVector(16) = "<<auxiliaryDerivativesVector(16)<<std::endl;
//    std::cout<<"auxiliaryDerivativesVector(31) = "<<auxiliaryDerivativesVector(31)<<std::endl;

//    std::cout<<"pow(auxiliaryEquationsVector(31),i-1) = "<<pow(auxiliaryEquationsVector(31),i-1)<<std::endl;
//    std::cout<<"temperaturePolyCoefficients(sectionT,i) = "<<temperaturePolyCoefficients(sectionT,i)<<std::endl;
//    std::cout<<"auxiliaryDerivativesVector(31)*i*pow(auxiliaryEquationsVector(31),i-1)*temperaturePolyCoefficients(sectionT,i) = "<<auxiliaryDerivativesVector(31)*i*pow(auxiliaryEquationsVector(31),i-1)*temperaturePolyCoefficients(sectionT,i)<<std::endl;
//    std::cout<<"This temperature loop is computed and the current temperature change = "<<auxiliaryDerivativesVector(34)<<std::endl;
//    std::cout<<" "<<std::endl;

};


    auxiliaryDerivativesVector(28) = auxiliaryDerivativesVector(30)*exp(auxiliaryEquationsVector(30)) ;                // u28

    if (referenceArea == 0.0){ // In case the density equation does weird things
        auxiliaryDerivativesVector(28) = 0.0;
    }


    auxiliaryDerivativesVector(33) = ((adiabeticIndex*specificGasConstant)/(2.0*auxiliaryEquationsVector(33)))*auxiliaryDerivativesVector(34);                // u33



    // Determine which section of the drag coefficient curve needs to be used

    for (int i=0; i < 5+1; i++){

        if (dragCoefficientMachRanges(i,0) <= auxiliaryEquationsVector(32) && auxiliaryEquationsVector(32) < dragCoefficientMachRanges(i,1)){

            sectionCD = i;
        }

    };

    auxiliaryDerivativesVector(32) = (auxiliaryEquationsVector(33)*auxiliaryDerivativesVector(15)-auxiliaryEquationsVector(15)*auxiliaryDerivativesVector(33))/(auxiliaryEquationsVector(33)*auxiliaryEquationsVector(33)); // u32

    auxiliaryDerivativesVector(29) = dragCoefficientPolyCoefficients(sectionCD,1)*auxiliaryDerivativesVector(32);              // u29




    auxiliaryDerivativesVector(27) = 0.5*referenceArea*auxiliaryEquationsVector(15)*(auxiliaryEquationsVector(15)*(auxiliaryEquationsVector(29)*auxiliaryDerivativesVector(28)+auxiliaryEquationsVector(28)*auxiliaryDerivativesVector(29))+
                                                                                     2.0*auxiliaryEquationsVector(28)*auxiliaryEquationsVector(29)*auxiliaryDerivativesVector(15));                // u27

//    // Accounting for altitude higher than 320 km MOLA: assumption, density can be neglected
//    if (auxiliaryEquationsVector(31)>=320.0){
//        auxiliaryDerivativesVector(27) = 0.0;
//        std::cout<<"It's done, it's gone!"<<std::endl;
//    }


// auxiliaryDerivativesVector() = ;                // u

    // Set vertical ascent to false again
    //        verticalInertialFlightPathAngleSet = false;
//            verticalInertialFlightPathAngleSet = NULL;
//            std::cout<<"verticalInertialFlightPathAngleSet der 3 = "<<verticalInertialFlightPathAngleSet<<std::endl;
    //        verticalRotationalFlightPathAngleSet = false;
//            verticalRotationalFlightPathAngleSet = NULL;
//            std::cout<<"verticalRotationalFlightPathAngleSet der 3 = "<<verticalRotationalFlightPathAngleSet<<std::endl;


    return auxiliaryDerivativesVector;
}

//////////////////////////////////////////////// Auxiliary Functions //////////////////////////////////////////////////////////////////////

Eigen::MatrixXd getAuxiliaryFunctions( const tudat::basic_mathematics::Vector7d& aState, const double time, const Eigen::Vector3d& thrustAccelerationsBframe, const Eigen::VectorXd& auxiliaryEquationsVector,
                                         const Eigen::VectorXd& auxiliaryDerivativesVector){

    auxiliaryFunctionsMatrix = Eigen::MatrixXd::Zero(35,39);       // Setting the complete matrix and filling it with zeros for now

//    Eigen::MatrixXd thrustAzimuthMatrix = MAV.thrustAzimuth();      // Setting the thrust angle matrices
//    Eigen::MatrixXd thrustElevationMatrix = MAV.thrustElevation();

    // The following expressions are described in the order in which the equations have to be computed corresponding to the respective vector entry
    // Which in this case means that the not all positions in the matrix will be used. The other values will simply be 0.

    // Avoid cosine rounding errors
//    double cx10x11;

//    if (abs(cos(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11)))<6.2e-17){
//        cx10x11 = 0;
//    }
//    else {
//        cx10x11 = cos(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11));
//    }


    // Determining the Thrust azimuth (psiT) and the Thrust elevation (epsilonT) angles based on the altitude
            // Determine the proper azimuth value for the current altitude section
            int sectionThrustAz = 0;    // Set the current azimuth value to the default first section
            for (int i = 0; i < thrustAzimuthMatrix.rows();i++){
                if (thrustAzimuthMatrix(i,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < thrustAzimuthMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                    sectionThrustAz = i;
//                    std::cout<<"sectionThrustAz = "<<sectionThrustAz+1<<std::endl;
                }
            }

                const double thrustAzimuth = thrustAzimuthMatrix(sectionThrustAz,2); // Set the thrust azimuth to the current azimuth corresponding to the current altitude section

            // Determine the proper elevation value for the current altitude section
            int sectionThrustEl = 0;    // Set the current elevation value to the default first section
            for (int i = 0; i < thrustElevationMatrix.rows();i++){
                if (thrustElevationMatrix(i,0) <= auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31) < thrustElevationMatrix(i,1)){ // Test for all the sections (independent of how many sections there are)
                    sectionThrustEl = i;
                }
            }

                const double thrustElevation = thrustElevationMatrix(sectionThrustEl,2); // Set the thrust elevation to the current elevation corresponding to the current altitude section

//            const double thrustAzimuth = 0.0; // Test
//            const double thrustElevation = 0.0; // Test


    double cx12;

    if (abs(cos(auxiliaryEquationsVector(12)))<6.2e-17){
        cx12 = 0;
    }
    else {
        cx12 = cos(auxiliaryEquationsVector(12));
    }

    double cx13;

    if (abs(cos(auxiliaryEquationsVector(13)))<6.2e-17){
        cx13 = 0;
    }
    else {
        cx13 = cos(auxiliaryEquationsVector(13));
    }

    double cx14;

    if (abs(cos(auxiliaryEquationsVector(14)))<6.2e-17){
        cx14 = 0;

    }
    else {
        cx14 = cos(auxiliaryEquationsVector(14));
    }

    // Same for sin

//   double sx10x11;

//    if (abs(sin(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11)))<1.2e-16){ // Using x10 and x11

//        sx10x11 = 0;

//    }
//    else {
//        sx10x11 = sin(auxiliaryEquationsVector(10)+auxiliaryEquationsVector(11));
//    }

    double sx12;

    if (abs(sin(auxiliaryEquationsVector(12)))<6.2e-17){
        sx12 = 0;
    }
    else {
        sx12 = sin(auxiliaryEquationsVector(12));
    }

    double sx13;

    if (abs(sin(auxiliaryEquationsVector(13)))<6.2e-17){
        sx13 = 0;
    }
    else {
        sx13 = sin(auxiliaryEquationsVector(13));
    }

    double sx14;

    if (abs(sin(auxiliaryEquationsVector(14)))<6.2e-17){
        sx14 = 0;
//        std::cout<<"cx14 has been set to 0"<<std::endl;
    }
    else {
        sx14 = sin(auxiliaryEquationsVector(14));
    }

    // w4
//    auxiliaryFunctionsMatrix(4,0) = auxiliaryEquationsVector(27)/auxiliaryEquationsVector(7);   // Added because of the mistake found in the recurrence relation of W4,2
//    auxiliaryFunctionsMatrix(4,1) = auxiliaryEquationsVector(1)/auxiliaryEquationsVector(9);
//    auxiliaryFunctionsMatrix(4,2) = auxiliaryEquationsVector(0);
//    auxiliaryFunctionsMatrix(4,3) = cx10x11;
//    auxiliaryFunctionsMatrix(4,3) = cos(auxiliaryEquationsVector(49));
    // Avoid cosine rounding errors
//            if (abs(auxiliaryFunctionsMatrix(4,3))<6.2e-17){
//              auxiliaryFunctionsMatrix(4,3) = 0;
//            }
    auxiliaryFunctionsMatrix(4,4) = sx12;
//    std::cout<<"w4,4 (sin(x12)) = "<<auxiliaryFunctionsMatrix(4,4)<<std::endl;
    auxiliaryFunctionsMatrix(4,5) = cx13;
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,5))<6.2e-17){
              auxiliaryFunctionsMatrix(4,5) = 0;
            }
    auxiliaryFunctionsMatrix(4,6) = cx12;
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,6))<6.2e-17){
              auxiliaryFunctionsMatrix(4,6) = 0;
            }
//    std::cout<<"w4,6 (cos(x12)) = "<<auxiliaryFunctionsMatrix(4,6)<<std::endl;
    auxiliaryFunctionsMatrix(4,7) = sx14;

    auxiliaryFunctionsMatrix(4,38) = cx14;
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,38))<6.2e-17){
              auxiliaryFunctionsMatrix(4,38) = 0;
            }
//    auxiliaryFunctionsMatrix(4,8) = sx10x11;
//    auxiliaryFunctionsMatrix(4,8) = sin(auxiliaryEquationsVector(49));
    // Avoid sine rounding errors
//    if (abs(auxiliaryFunctionsMatrix(4,8))<1.22e-16){
//        auxiliaryFunctionsMatrix(4,8) = 0;
//    }
    /// Debug ///

//    std::cout<<"sin(x13) = "<<sin(auxiliaryEquationsVector(13))<<std::endl;
//    std::cout<<"sx13 = "<<sx13<<std::endl;

    /// Debug ///
    auxiliaryFunctionsMatrix(4,9) = sx13;
    auxiliaryFunctionsMatrix(4,10) = auxiliaryFunctionsMatrix(4,4)*auxiliaryFunctionsMatrix(4,5);
    auxiliaryFunctionsMatrix(4,11) = auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,7);
//    auxiliaryFunctionsMatrix(4,12) = auxiliaryFunctionsMatrix(4,9)*auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(4,12) = auxiliaryFunctionsMatrix(4,9)*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(4,13) = auxiliaryFunctionsMatrix(4,4)*auxiliaryFunctionsMatrix(4,9);
    auxiliaryFunctionsMatrix(4,14) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,5);
//    auxiliaryFunctionsMatrix(4,15) = auxiliaryFunctionsMatrix(4,6)*auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(4,15) = auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(4,16) = auxiliaryFunctionsMatrix(4,9)*auxiliaryFunctionsMatrix(4,7);
//    auxiliaryFunctionsMatrix(4,17) = auxiliaryFunctionsMatrix(4,10)*auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(4,17) = auxiliaryFunctionsMatrix(4,10)*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(4,18) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,12);
    auxiliaryFunctionsMatrix(4,19) = auxiliaryFunctionsMatrix(4,3)*auxiliaryFunctionsMatrix(4,13);
    auxiliaryFunctionsMatrix(4,20) = auxiliaryFunctionsMatrix(4,10)*auxiliaryFunctionsMatrix(4,7);
    auxiliaryFunctionsMatrix(4,21) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,16);
    auxiliaryFunctionsMatrix(4,22) = auxiliaryFunctionsMatrix(4,3)*(auxiliaryFunctionsMatrix(4,11)-auxiliaryFunctionsMatrix(4,17));
    auxiliaryFunctionsMatrix(4,23) = auxiliaryFunctionsMatrix(4,3)*(-auxiliaryFunctionsMatrix(4,20)-auxiliaryFunctionsMatrix(4,15));
    auxiliaryFunctionsMatrix(4,24) = auxiliaryFunctionsMatrix(4,2)*(auxiliaryFunctionsMatrix(4,22)-auxiliaryFunctionsMatrix(4,18));

    auxiliaryFunctionsMatrix(4,25) = Thrust/auxiliaryEquationsVector(7);
    auxiliaryFunctionsMatrix(4,26) = cos(thrustAzimuth);
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,26))<6.2e-17){
              auxiliaryFunctionsMatrix(4,26) = 0;
            }

    auxiliaryFunctionsMatrix(4,27) = cos(thrustElevation);
    // Avoid cosine rounding errors
            if (abs(auxiliaryFunctionsMatrix(4,27))<6.2e-17){
              auxiliaryFunctionsMatrix(4,27) = 0;
            }

    auxiliaryFunctionsMatrix(4,28) = sin(thrustAzimuth);
    auxiliaryFunctionsMatrix(4,29) = sin(thrustElevation);
    auxiliaryFunctionsMatrix(4,30) = auxiliaryFunctionsMatrix(4,26)*auxiliaryFunctionsMatrix(4,27);
    auxiliaryFunctionsMatrix(4,31) = auxiliaryFunctionsMatrix(4,28)*auxiliaryFunctionsMatrix(4,27);
    auxiliaryFunctionsMatrix(4,32) = auxiliaryFunctionsMatrix(4,25)*auxiliaryFunctionsMatrix(4,30);
    auxiliaryFunctionsMatrix(4,33) = auxiliaryFunctionsMatrix(4,25)*auxiliaryFunctionsMatrix(4,31);
    auxiliaryFunctionsMatrix(4,34) = auxiliaryFunctionsMatrix(4,25)*auxiliaryFunctionsMatrix(4,29);
    auxiliaryFunctionsMatrix(4,35) = auxiliaryFunctionsMatrix(4,33)*(auxiliaryFunctionsMatrix(4,19)-auxiliaryFunctionsMatrix(4,14));
    auxiliaryFunctionsMatrix(4,36) = auxiliaryFunctionsMatrix(4,34)*(auxiliaryFunctionsMatrix(4,23)-auxiliaryFunctionsMatrix(4,21));
    auxiliaryFunctionsMatrix(4,37) = 1.0/auxiliaryEquationsVector(7);



    // w5
//    auxiliaryFunctionsMatrix(5,1) = auxiliaryEquationsVector(2)/auxiliaryEquationsVector(9);
    auxiliaryFunctionsMatrix(5,2) = auxiliaryFunctionsMatrix(4,8)*(auxiliaryFunctionsMatrix(4,11)-auxiliaryFunctionsMatrix(4,17));
    auxiliaryFunctionsMatrix(5,3) = auxiliaryFunctionsMatrix(4,3)*auxiliaryFunctionsMatrix(4,12);
    auxiliaryFunctionsMatrix(5,4) = auxiliaryFunctionsMatrix(4,8)*auxiliaryFunctionsMatrix(4,13);
    auxiliaryFunctionsMatrix(5,5) = auxiliaryFunctionsMatrix(4,3)*auxiliaryFunctionsMatrix(4,5);
    auxiliaryFunctionsMatrix(5,6) = auxiliaryFunctionsMatrix(4,8)*(-auxiliaryFunctionsMatrix(4,20)-auxiliaryFunctionsMatrix(4,11));
    auxiliaryFunctionsMatrix(5,7) = auxiliaryFunctionsMatrix(4,3)*auxiliaryFunctionsMatrix(4,16);
    auxiliaryFunctionsMatrix(5,8) = auxiliaryFunctionsMatrix(4,2)*(auxiliaryFunctionsMatrix(5,2)+auxiliaryFunctionsMatrix(5,3));
    auxiliaryFunctionsMatrix(5,9) = auxiliaryFunctionsMatrix(4,33)*(auxiliaryFunctionsMatrix(5,4)+auxiliaryFunctionsMatrix(5,5));
    auxiliaryFunctionsMatrix(5,10) = auxiliaryFunctionsMatrix(4,34)*(auxiliaryFunctionsMatrix(5,6)+auxiliaryFunctionsMatrix(5,7));


    // w6
    auxiliaryFunctionsMatrix(6,0) = auxiliaryFunctionsMatrix(4,4)*auxiliaryFunctionsMatrix(4,7);  // Added because of the mistake found in the complete transformation matrix
//    auxiliaryFunctionsMatrix(6,1) = auxiliaryEquationsVector(3)/auxiliaryEquationsVector(9);
    auxiliaryFunctionsMatrix(6,2) = auxiliaryFunctionsMatrix(4,5)*auxiliaryFunctionsMatrix(4,15);
    auxiliaryFunctionsMatrix(6,3) = auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,9);
    auxiliaryFunctionsMatrix(6,4) = auxiliaryFunctionsMatrix(4,5)*auxiliaryFunctionsMatrix(4,11);
//    auxiliaryFunctionsMatrix(6,5) = auxiliaryFunctionsMatrix(4,4)*auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(6,5) = auxiliaryFunctionsMatrix(4,4)*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(6,6) = auxiliaryFunctionsMatrix(4,2)*(auxiliaryFunctionsMatrix(6,2)+auxiliaryFunctionsMatrix(6,0)); // Changed becuase of the mistake found in the complete transformation matrix
    auxiliaryFunctionsMatrix(6,7) = auxiliaryFunctionsMatrix(4,33)*auxiliaryFunctionsMatrix(6,3);
    auxiliaryFunctionsMatrix(6,8) = auxiliaryFunctionsMatrix(4,34)*(auxiliaryFunctionsMatrix(6,4)-auxiliaryFunctionsMatrix(6,5));



//    // w8
//    auxiliaryFunctionsMatrix(8,1) = auxiliaryEquationsVector(1)*auxiliaryEquationsVector(4);
//    auxiliaryFunctionsMatrix(8,2) = auxiliaryEquationsVector(2)*auxiliaryEquationsVector(5);
//    auxiliaryFunctionsMatrix(8,3) = auxiliaryEquationsVector(3)*auxiliaryEquationsVector(6);





    // w9
//    auxiliaryFunctionsMatrix(9,1) = (auxiliaryEquationsVector(9)*auxiliaryDerivativesVector(8))/auxiliaryEquationsVector(8);






    // w11
    auxiliaryFunctionsMatrix(11,0) = auxiliaryEquationsVector(15)*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(11,1) = auxiliaryFunctionsMatrix(11,0)/auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(11,2) = auxiliaryFunctionsMatrix(11,1)*auxiliaryFunctionsMatrix(4,9);
    // Avoid singularities
    if (auxiliaryFunctionsMatrix(4,6) == 0.0){
        auxiliaryFunctionsMatrix(11,3) = 0.0;
    }
    else {
    auxiliaryFunctionsMatrix(11,3) = auxiliaryFunctionsMatrix(11,2)/auxiliaryFunctionsMatrix(4,6);
    }

    // w12
    auxiliaryFunctionsMatrix(12,1) = auxiliaryFunctionsMatrix(11,1)*auxiliaryFunctionsMatrix(4,5);

    // w13
    auxiliaryFunctionsMatrix(13,0) = auxiliaryFunctionsMatrix(4,28)*auxiliaryFunctionsMatrix(4,27);
    auxiliaryFunctionsMatrix(13,1) = auxiliaryFunctionsMatrix(4,7)*auxiliaryFunctionsMatrix(4,5);
    auxiliaryFunctionsMatrix(13,2) = rotationalVelocity*rotationalVelocity*auxiliaryEquationsVector(16)*auxiliaryFunctionsMatrix(4,6);
    auxiliaryFunctionsMatrix(13,3) = auxiliaryFunctionsMatrix(13,2)*auxiliaryFunctionsMatrix(4,4);
    auxiliaryFunctionsMatrix(13,4) = Thrust*auxiliaryFunctionsMatrix(13,0)*auxiliaryFunctionsMatrix(4,37);
    auxiliaryFunctionsMatrix(13,5) = auxiliaryFunctionsMatrix(13,3)*auxiliaryFunctionsMatrix(4,9)+auxiliaryFunctionsMatrix(13,4);
    // Avoid singularities
    if (auxiliaryEquationsVector(15) == 0.0){
        auxiliaryFunctionsMatrix(13,6) = 0.0;
    }
    else {
    auxiliaryFunctionsMatrix(13,6) = auxiliaryFunctionsMatrix(13,5)/auxiliaryEquationsVector(15);
    }
    auxiliaryFunctionsMatrix(13,7) = -2.0*rotationalVelocity*auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(13,1)+auxiliaryFunctionsMatrix(13,6);
    // Avoid singularities
    if (auxiliaryFunctionsMatrix(4,38) == 0.0){
        auxiliaryFunctionsMatrix(13,8) = 0.0;
    }
    else {
        auxiliaryFunctionsMatrix(13,8) = auxiliaryFunctionsMatrix(13,7)/auxiliaryFunctionsMatrix(4,38);
    }
    auxiliaryFunctionsMatrix(13,9) = (2.0*rotationalVelocity+auxiliaryFunctionsMatrix(11,3))*auxiliaryFunctionsMatrix(4,4)+auxiliaryFunctionsMatrix(13,8);

    // w14
    auxiliaryFunctionsMatrix(14,0) = Thrust*auxiliaryFunctionsMatrix(4,29)*auxiliaryFunctionsMatrix(4,37);
    auxiliaryFunctionsMatrix(14,1) = auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,38)+auxiliaryFunctionsMatrix(13,1)*auxiliaryFunctionsMatrix(4,4);
    auxiliaryFunctionsMatrix(14,2) = -standardGravitationalParameter*auxiliaryFunctionsMatrix(4,38);
    auxiliaryFunctionsMatrix(14,3) = auxiliaryEquationsVector(16)*auxiliaryEquationsVector(16);
    auxiliaryFunctionsMatrix(14,4) = auxiliaryFunctionsMatrix(14,2)/auxiliaryFunctionsMatrix(14,3);
    auxiliaryFunctionsMatrix(14,5) = auxiliaryFunctionsMatrix(14,4)+auxiliaryFunctionsMatrix(13,2)*auxiliaryFunctionsMatrix(14,1)+auxiliaryFunctionsMatrix(14,0);
    // Avoid singularities
    if (auxiliaryEquationsVector(15) == 0.0){
        auxiliaryFunctionsMatrix(14,6) = 0.0;
    }
    else {
    auxiliaryFunctionsMatrix(14,6) = auxiliaryFunctionsMatrix(14,5)/auxiliaryEquationsVector(15);
    }

    auxiliaryFunctionsMatrix(14,7) = 2.0*rotationalVelocity*auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,9)+auxiliaryFunctionsMatrix(11,1)+auxiliaryFunctionsMatrix(14,6);

    /// Debug ///
//    std::cout<<"w14,5 = "<<auxiliaryFunctionsMatrix(14,5)<<std::endl;
//    std::cout<<"x15 = "<<auxiliaryEquationsVector(15)<<std::endl;
//    std::cout<<"w4,6 = "<<auxiliaryFunctionsMatrix(4,6)<<std::endl;
//    std::cout<<"w4,9 = "<<auxiliaryFunctionsMatrix(4,9)<<std::endl;
//    std::cout<<"2.0*rotationalVelocity*auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,9) = "<<2.0*rotationalVelocity*auxiliaryFunctionsMatrix(4,6)*auxiliaryFunctionsMatrix(4,9)<<std::endl;
//    std::cout<<"w11,1 = "<<auxiliaryFunctionsMatrix(11,1)<<std::endl;
//    std::cout<<"w14,6 = "<<auxiliaryFunctionsMatrix(14,6)<<std::endl;
//    std::cout<<"w14,7 = "<<auxiliaryFunctionsMatrix(14,7)<<std::endl;
    /// Debug ///


    // w15
    auxiliaryFunctionsMatrix(15,0) = Thrust*auxiliaryFunctionsMatrix(4,26)*auxiliaryFunctionsMatrix(4,27)-auxiliaryEquationsVector(27);
    auxiliaryFunctionsMatrix(15,1) = auxiliaryFunctionsMatrix(15,0)/auxiliaryEquationsVector(7);
    auxiliaryFunctionsMatrix(15,2) = -standardGravitationalParameter*auxiliaryFunctionsMatrix(4,7);
    auxiliaryFunctionsMatrix(15,3) = auxiliaryFunctionsMatrix(15,2)/auxiliaryFunctionsMatrix(14,3);
    auxiliaryFunctionsMatrix(15,4) = auxiliaryFunctionsMatrix(4,38)*auxiliaryFunctionsMatrix(4,5);
    auxiliaryFunctionsMatrix(15,5) = auxiliaryFunctionsMatrix(4,7)*auxiliaryFunctionsMatrix(4,6)-auxiliaryFunctionsMatrix(15,4)*auxiliaryFunctionsMatrix(4,4);
    auxiliaryFunctionsMatrix(15,6) = auxiliaryFunctionsMatrix(13,2)*auxiliaryFunctionsMatrix(15,5)+auxiliaryFunctionsMatrix(15,1)+auxiliaryFunctionsMatrix(15,3);

    /// Debug ///
//    std::cout<<"u15 = "<<auxiliaryFunctionsMatrix(15,6)<<std::endl;
    /// Debug ///

    // w16
    auxiliaryFunctionsMatrix(16,1) = auxiliaryEquationsVector(15)*auxiliaryFunctionsMatrix(4,7);



    // w27
    auxiliaryFunctionsMatrix(27,1) = auxiliaryEquationsVector(29)*auxiliaryDerivativesVector(28);
    auxiliaryFunctionsMatrix(27,2) = auxiliaryEquationsVector(28)*auxiliaryDerivativesVector(29);
    auxiliaryFunctionsMatrix(27,3) = auxiliaryEquationsVector(28)*auxiliaryEquationsVector(29);
    auxiliaryFunctionsMatrix(27,4) = auxiliaryEquationsVector(15)*(auxiliaryFunctionsMatrix(27,1)+auxiliaryFunctionsMatrix(27,2));
    auxiliaryFunctionsMatrix(27,5) = auxiliaryFunctionsMatrix(27,3)*auxiliaryDerivativesVector(15);
    auxiliaryFunctionsMatrix(27,6) = auxiliaryEquationsVector(15)*(auxiliaryFunctionsMatrix(27,4)+auxiliaryFunctionsMatrix(27,5));

//    // Accounting for altitude higher than 320 km MOLA: assumption, density can be neglected
//    if (auxiliaryEquationsVector(31)>=320.0){
//        auxiliaryDerivativesVector(27) = 0.0;
//        std::cout<<"It's done, it's gone!"<<std::endl;
//    }







    // w28
    auxiliaryFunctionsMatrix(28,1) = auxiliaryDerivativesVector(30)*auxiliaryEquationsVector(28);






    // w30
    auxiliaryFunctionsMatrix(30,1) = pow(auxiliaryEquationsVector(31),9.0);
    auxiliaryFunctionsMatrix(30,2) = pow(auxiliaryEquationsVector(31),8.0);
    auxiliaryFunctionsMatrix(30,3) = pow(auxiliaryEquationsVector(31),7.0);
    auxiliaryFunctionsMatrix(30,4) = pow(auxiliaryEquationsVector(31),6.0);
    auxiliaryFunctionsMatrix(30,5) = pow(auxiliaryEquationsVector(31),5.0);
    auxiliaryFunctionsMatrix(30,6) = pow(auxiliaryEquationsVector(31),4.0);
    auxiliaryFunctionsMatrix(30,7) = pow(auxiliaryEquationsVector(31),3.0);
    auxiliaryFunctionsMatrix(30,8) = pow(auxiliaryEquationsVector(31),2.0);

    for (int i=2; i<10+1; i++){

        if (i==2){

//            auxiliaryFunctionsMatrix(30,9) = auxiliaryDerivativesVector(31)*(2.0*densityPolyCoefficients(2)*auxiliaryEquationsVector(31)+densityPolyCoefficients(1));
            auxiliaryFunctionsMatrix(30,19) = (2.0*densityPolyCoefficients(2)*auxiliaryEquationsVector(31));
        }
        else {

//           auxiliaryFunctionsMatrix(30,9) += auxiliaryDerivativesVector(31)*i*densityPolyCoefficients(i)*auxiliaryFunctionsMatrix(30,(11-i));
            auxiliaryFunctionsMatrix(30,19) += i*densityPolyCoefficients(i)*auxiliaryFunctionsMatrix(30,(11-i));

        };
    };
    auxiliaryFunctionsMatrix(30,9) = auxiliaryDerivativesVector(31)*auxiliaryFunctionsMatrix(30,19)+densityPolyCoefficients(1)*auxiliaryFunctionsMatrix(31);






    // w32
    auxiliaryFunctionsMatrix(32,1) = auxiliaryEquationsVector(33)*auxiliaryDerivativesVector(15);
    auxiliaryFunctionsMatrix(32,2) = auxiliaryEquationsVector(15)*auxiliaryDerivativesVector(33);
    auxiliaryFunctionsMatrix(32,3) = auxiliaryEquationsVector(33)*auxiliaryEquationsVector(33);
    auxiliaryFunctionsMatrix(32,4) = (auxiliaryFunctionsMatrix(32,1)-auxiliaryFunctionsMatrix(32,2))/(auxiliaryFunctionsMatrix(32,3));





    // w33
    auxiliaryFunctionsMatrix(33,1) = auxiliaryDerivativesVector(34)/auxiliaryEquationsVector(33);




    // w34
    // First it has to be determined whether or not these functions have to be computed. If not, they remain zero. If so they only one of them is computed.

    if (temperatureAltitudeRanges(1,0)<=auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31)<temperatureAltitudeRanges(1,1)){   // Section 2

        auxiliaryFunctionsMatrix(34,12) = (2.0*temperaturePolyCoefficients(1,2)*auxiliaryEquationsVector(31));

      auxiliaryFunctionsMatrix(34,2) = auxiliaryDerivativesVector(31)*auxiliaryFunctionsMatrix(34,12)+temperaturePolyCoefficients(1,1)*auxiliaryDerivativesVector(31);
    }
    else if (temperatureAltitudeRanges(2,0)<=auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31)<temperatureAltitudeRanges(2,1)){ // Section 3

        for (int i=2; i<6+1;i++){

            if (i==2){

//                auxiliaryFunctionsMatrix(34,3) = auxiliaryDerivativesVector(31)*(2.0*temperaturePolyCoefficients(2,2)*auxiliaryEquationsVector(31)+temperaturePolyCoefficients(2,1));
                auxiliaryFunctionsMatrix(34,13) = (2.0*temperaturePolyCoefficients(2,2)*auxiliaryEquationsVector(31));
            }


            else {

//                auxiliaryFunctionsMatrix(34,3) += auxiliaryDerivativesVector(31)*i*temperaturePolyCoefficients(2,i)*auxiliaryFunctionsMatrix(30,(11-i));
                auxiliaryFunctionsMatrix(34,13) += i*temperaturePolyCoefficients(2,i)*auxiliaryFunctionsMatrix(30,(11-i));
//                std::cout<<"i*temperaturePolyCoefficients(2,i)*auxiliaryFunctionsMatrix(30,(11-i)) = "<<i*temperaturePolyCoefficients(2,i)*auxiliaryFunctionsMatrix(30,(11-i))<<std::endl;
            };


        };

    auxiliaryFunctionsMatrix(34,3) = auxiliaryDerivativesVector(31)*auxiliaryFunctionsMatrix(34,13)+temperaturePolyCoefficients(2,1)*auxiliaryDerivativesVector(31);
//    std::cout<<"auxiliaryFunctionsMatrix(34,3) = "<<auxiliaryFunctionsMatrix(34,3)<<std::endl;

    }
    else if (temperatureAltitudeRanges(3,0)<=auxiliaryEquationsVector(31) && auxiliaryEquationsVector(31)<temperatureAltitudeRanges(3,1)){ // Section 4

        for (int i=2; i<8+1;i++){

            if (i==2){

//                auxiliaryFunctionsMatrix(34,4) = auxiliaryDerivativesVector(31)*(2.0*temperaturePolyCoefficients(3,2)*auxiliaryEquationsVector(31)+temperaturePolyCoefficients(3,1));
                auxiliaryFunctionsMatrix(34,14) = (2.0*temperaturePolyCoefficients(3,2)*auxiliaryEquationsVector(31));
            }


            else {

//                auxiliaryFunctionsMatrix(34,4) += auxiliaryDerivativesVector(31)*i*temperaturePolyCoefficients(3,i)*auxiliaryFunctionsMatrix(30,(11-i));
                auxiliaryFunctionsMatrix(34,14) += i*temperaturePolyCoefficients(3,i)*auxiliaryFunctionsMatrix(30,(11-i));

            };

        };


        auxiliaryFunctionsMatrix(34,4) = auxiliaryDerivativesVector(31)*auxiliaryFunctionsMatrix(34,14)+temperaturePolyCoefficients(3,1)*auxiliaryDerivativesVector(31);
};




// auxiliaryFunctionsMatrix(,)


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
 Eigen::MatrixXd thrustAzimuthPolyCoefficients;                 // P_psiT   these are the polynomial coefficients for the fit for the psiT curve
 Eigen::MatrixXd thrustElevationPolyCoefficients;               // P_epsilonT   these are the polynomial coefficients for the fit for the epsilonT curve

    // Set functions

  double FlightPathAngle;         // Flight path angle in rad
 double HeadingAngle;            // Heading angle in rad

// double rotationalFlightPathAngle;         // Rotational flight path angle in rad
// double inertialFlightPathAngle;           // Inertial flight path angle in rad
// double rotationalHeadingAngle;            // Rotational heading angle in rad
// double inertialHeadingAngle;              // Inertial heading angle in rad

// bool rotationalFlightPathAngleSet;         // All of these are used to let the program know that a predefined angle was set and that that angle should be used (initially)
// bool inertialFlightPathAngleSet;
// bool rotationalHeadingAngleSet;
// bool inertialHeadingAngleSet;


// bool verticalRotationalFlightPathAngleSet;       // All of these are used for the vertical ascent case
// bool verticalInertialFlightPathAngleSet;
// bool verticalRotationalHeadingAngleSet;
// bool verticalInertialHeadingAngleSet;





    // Additional in-class used variables


    int sectionT;                                   // This variable holds the "(section number -1)" for the temperature curve fit
    int powerT;                                     // This variable holds the section corresponding order for the temperature curve fit

    int sectionCD;                                  // This variable holds the "(section number -1)" for the drag coefficient curve fit

    Eigen::VectorXd auxiliaryEquationsVector;               // The vector containing the auxiliary equations xn
    Eigen::VectorXd auxiliaryDerivativesVector;             // The vector containing the auxiliary derivatives un
    Eigen::MatrixXd auxiliaryFunctionsMatrix;               // The matrix containing the auxiliary functions wn,m

};

#endif // AUXILIARY_H
