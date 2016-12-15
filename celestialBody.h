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
 *      160630    S.D. Petrovic     Found mistake in molar mass of Martian atmosphere (typo) and corrected it
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H


#include <iostream>
#include <Eigen/Core>
#include <string>


// Author: Stacha Petrovic
// Date created: 11th April 2016


class celestialBody

        /* This class will describe the different celestial body characteristics
         * including atmospheric data.
         * These constants are:
         *
         * - gamma_a = 1.35                                         adiabetic index for Mars
         * - Rstar = Ra/Ma =~ 191.842635            [m^2/(s^2*K)]   specific gas constant for Mars (Ra = 8.3144598 m^2*kg/(s^2*K*mol) and Ma = 0.04334 kg/mol)
         * - PTn                                    [-]             these are the polynomial coefficients for the fit for the temperature curves
         * - Prho n                                 [-]             these are the polynomial coefficients for the fit for the density curves
         * - mu_M = 4.283*10^13                     [m^3/s^2]       standard gravitational parameter for Mars
         * - marsRotationalVelocity = 7.088*10^-5   [rad/s]         rotational velocity of Mars
         * - OmegaP = 0                             [rad]           relative angle between the prime meridian and the x-axis (for now it is chosen such that
         *                                                          the inertial frame was set at time t = start time and going through the prime meridian)
         * - t0 = 0                                 [s]             time between the start time and the time that the inertial frame was set
         * - temperatureAltitudeRanges              [km MOLA]       altitude range per section for the temperature-altitude curve
         * - Rm = 3396e3                              [m]             MOLA radius of Mars
         *
         */





{
public:

    // Setting the celestial body with the default being Mars (is the constructor)

    celestialBody(const std::string planet = "Mars"){

        std::cout<<"You have chosen "<<planet<<" as your celestial body"<<std::endl;

        if (planet == "Mars" || planet == "mars")
        {

            /// Set different Mars constants: ///


            // Adiabetic Index
                adiabeticIndex_ = 1.35;


            // Specific Gas Constant
                const double massOfMars = 0.04334;         // [kg/mol]
//                const double gasConstant = 8.3144598;       // [m^2*kg/(s^2*K*mol)]
                const double gasConstant = 8.3144598e-6;    // [km^2*kg/(s^2*K*mol)]

//                specificGasConstant_ = gasConstant/massOfMars;  // [m^2/(s^2*K)]
                specificGasConstant_ = gasConstant/massOfMars;  // [km^2/(s^2*K)]


            // Standard gravitational parameter
//                standardGravitationalParameter_ = 4.2828314e13;  // [m^3/s^2]
                standardGravitationalParameter_ = 4.2828314e4;  // [km^3/s^2]


            // Rotational velocity
               rotationalVelocity_ = 7.088e-5;  // [rad/s]
//                rotationalVelocity_ = 0;// [rad/s] for verification non-rotating planet...


            // Prime Meridian angle
                primeMeridianAngle_ = 0.0;  // [rad]  (please note that this angle is defined negative in the positive direction when looking from the inertial frame to the rotating frame. I might change this later...)


            // Inertial Frame Time
                inertialFrameTime_ = 0.0;  // [s]

            // MOLA radius of Mars
//                bodyReferenceRadius_ = 3396e3; // [m]
                bodyReferenceRadius_ = 3396.0; // [km]


            // Temperature Polynomial Coefficients
                temperaturePolyCoefficients_ = Eigen::MatrixXd::Zero(5, 9);  // PTn      Temperature polynomial coefficients


                // Section 1: -0.6 to 5.04 km MOLA
                temperaturePolyCoefficients_(0,0) = 1.941652034325321e2; //194.165;
                temperaturePolyCoefficients_(0,1) = 3.415126502921027; //3.415;

                // Section 2: 5.04 to 35.53 km MOLA
                temperaturePolyCoefficients_(1,0) = 2.220519114287426e2; //222.051;
                temperaturePolyCoefficients_(1,1) = -2.129543575042308; //-2.130;
                temperaturePolyCoefficients_(1,2) = 0.005776725906697; //0.006;

                // Section 3: 35.53 to 75.07 km MOLA
                temperaturePolyCoefficients_(2,0) = -1.166700024356974e4; //-1.167e4;
                temperaturePolyCoefficients_(2,1) = 1.407033824881397e3; //1.407e3;
                temperaturePolyCoefficients_(2,2) = -68.293490289773885; //-68.294;
                temperaturePolyCoefficients_(2,3) = 1.732604336852450; //1.733;
                temperaturePolyCoefficients_(2,4) = -0.024282109172326; //-0.0243;
                temperaturePolyCoefficients_(2,5) = 1.785414835183911e-4; //1.785e-4;
                temperaturePolyCoefficients_(2,6) = -5.388211105876920e-7; //-5.388e-7;

                // Section 4: 75.07 to 170.05 km MOLA
                temperaturePolyCoefficients_(3,0) = 2.235938205082397e5; //2.236e5;
                temperaturePolyCoefficients_(3,1) = -1.522453634203182e4; //-1.523e4;
                temperaturePolyCoefficients_(3,2) = 4.473778178802255e2; //447.378;
                temperaturePolyCoefficients_(3,3) = -7.404992496139436; //-7.405;
                temperaturePolyCoefficients_(3,4) = 0.075523662530184; //0.076;
                temperaturePolyCoefficients_(3,5) = -4.862260819979541e-4; //-4.862e-4;
                temperaturePolyCoefficients_(3,6) = 1.930903736685301e-6; //1.931e-6;
                temperaturePolyCoefficients_(3,7) = -4.327597099336157e-9; //-4.328e-9;
                temperaturePolyCoefficients_(3,8) = 4.194150230099975e-12; //4.1942e-12;

                // Section 5: 170.05 to 320.0 km MOLA
                temperaturePolyCoefficients_(4,0) = 136.5;


            // Altitude range per section for the temperature-altitude curve
                temperatureAltitudeRanges_ = Eigen::MatrixXd::Zero(5,2);

                // Section 1
                temperatureAltitudeRanges_(0,0) = -0.6;   //Lower bound
                temperatureAltitudeRanges_(0,1) = 5.04;   //Upper bound

                // Section 2
                temperatureAltitudeRanges_(1,0) = 5.04;   //Lower bound
                temperatureAltitudeRanges_(1,1) = 35.53;   //Upper bound

                // Section 3
                temperatureAltitudeRanges_(2,0) = 35.53;   //Lower bound
                temperatureAltitudeRanges_(2,1) = 75.07;   //Upper bound

                // Section 4
                temperatureAltitudeRanges_(3,0) = 75.07;   //Lower bound
                temperatureAltitudeRanges_(3,1) = 170.05;   //Upper bound

                // Section 5
                temperatureAltitudeRanges_(4,0) = 170.05;   //Lower bound
                temperatureAltitudeRanges_(4,1) = 320.0;   //Upper bound


            // Density Polynomial Coefficients
                densityPolyCoefficients_ = Eigen::VectorXd::Zero(11);  // Prho n Density polynomial coefficients

                densityPolyCoefficients_(0) = -4.172119763435839; //-4.172;
                densityPolyCoefficients_(1) = -0.096179484326282; //-0.0962;
                densityPolyCoefficients_(2) = 0.001413771036602; //1.414e-3;
                densityPolyCoefficients_(3) = -9.603504536940288e-5; //-9.604e-5;
                densityPolyCoefficients_(4) = 2.272923177465556e-6; //2.273e-6;
                densityPolyCoefficients_(5) = -2.883775814622154e-8; //-2.884e-8;
                densityPolyCoefficients_(6) = 2.145912177247706e-10; //2.146e-10;
                densityPolyCoefficients_(7) = -9.620258504409526e-13; //-9.620e-13;
                densityPolyCoefficients_(8) = 2.558809216119519e-15; //2.559e-15;
                densityPolyCoefficients_(9) = -3.724118893425651e-18; //-3.724e-18;
                densityPolyCoefficients_(10) = 2.287447109373264e-21; //2.287e-21;





        }

        else{

            // Setting all paramters to 0

            adiabeticIndex_= 0.0;                                         // gamma_a      adiabetic index
            specificGasConstant_ = 0.0;                                   // Rstar        specific gas constant
            standardGravitationalParameter_ = 0.0;                        // mu_M         standard gravitational parameter
            rotationalVelocity_ = 0.0;                                    // rotational velocity of celestial body
            primeMeridianAngle_ = 0.0;                                    // OmegaP       relative angle between the prime meridian and the x-axis
            inertialFrameTime_ = 0.0;                                     // t0           time between the start time and the time that the inertial frame was set
            bodyReferenceRadius_ = 0.0;                                        // Rm           MOLA radius of Mars
            temperaturePolyCoefficients_ = Eigen::MatrixXd::Zero(1,1);  // PTn    temperature polynomial coefficients
            temperatureAltitudeRanges_ = Eigen::MatrixXd::Zero(1,1);    // altitude range per section for the temperature-altitude curve
            densityPolyCoefficients_ = Eigen::VectorXd::Zero(1);        // Prho n density polynomial coefficients




            std::cerr<<planet<<" is either not a valid string or has not yet been defined"<<std::endl;
        };




                                                    } // End of constructor

    // Destructor:

    ~celestialBody(void){
//        std::cout<<"The planet has been deleted!"<<std::endl;

    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////// In order to add extra celestial bodies you simple add another elseif statement and provide the required information //////
    ///////// PLEASE NOTE THAT IF NOT ALL INFORMATION IS PROVIDED IT COULD BE THAT A VALUE CLOSE TO 0 IS PROVIDED INSTEAD!//////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set functions

    void setStandardGravitationalParameter(const double updatedStandardGravitationalParameter)  // This function can be used to set the value of the Standard gravitational parameter of the planet. If set to 0, gravity is not taken into account in the equations
    {
        standardGravitationalParameter_ = updatedStandardGravitationalParameter;
    }

    void setRotationalVelocity(const double rotationalVelocity)         // This function can be used to set the value of the rotational Velocity of the celestial body. For instance for verification purposes.
    {
        rotationalVelocity_ = rotationalVelocity;
    }

    void setUpperAltitudeBound(const double upperAltitudeBound)         // This function can be used to adjust the final altitude
    {
        temperatureAltitudeRanges_(4,1) = upperAltitudeBound;
    }


    // Returning the different constant parameters

    const double adiabeticIndex() { return adiabeticIndex_; }                                   // gamma_a      adiabetic index
    const double specificGasConstant() { return specificGasConstant_; }                         // Rstar        specific gas constant
    const double standardGravitationalParameter() { return standardGravitationalParameter_; }   // mu_M         standard gravitational parameter
    const double rotationalVelocity() { return rotationalVelocity_; }                           // rotational velocity of Mars
    const double primeMeridianAngle() { return primeMeridianAngle_; }                           // OmegaP       relative angle between the prime meridian and the x-axis
    const double inertialFrameTime() { return inertialFrameTime_; }                             // t0           time between the start time and the time that the inertial frame was set
    const double bodyReferenceRadius() { return bodyReferenceRadius_; }                                   // Rm           MOLA radius of Mars

    // Returning the different polynomial coefficient parameter matrices

    const Eigen::MatrixXd temperaturePolyCoefficients() { return temperaturePolyCoefficients_;} // PTn    temperature polynomial coefficients

    const Eigen::MatrixXd temperatureAltitudeRanges() { return temperatureAltitudeRanges_;}     // altitude range per section for the temperature-altitude curve

    const Eigen::VectorXd densityPolyCoefficients() { return densityPolyCoefficients_;}         // Prho n density polynomial coefficients











protected:

private:

    // Creating the different constant parameters

    double adiabeticIndex_;                         // gamma_a      adiabetic index
    double specificGasConstant_;                    // Rstar        specific gas constant
    double standardGravitationalParameter_;         // mu_M         standard gravitational parameter
    double rotationalVelocity_;                     // rotational velocity of celestial body
    double primeMeridianAngle_;                     // OmegaP       relative angle between the prime meridian and the x-axis
    double inertialFrameTime_;                      // t0           time between the start time and the time that the inertial frame was set
    double bodyReferenceRadius_;                         // Rm           MOLA radius of Mars

    // Creating the different polynomial coefficient parameter matrices

    Eigen::MatrixXd temperaturePolyCoefficients_;   // PTn    temperature polynomial coefficients

    Eigen::MatrixXd temperatureAltitudeRanges_;     // altitude range per section for the temperature-altitude curve

    Eigen::VectorXd densityPolyCoefficients_;       // Prho n density polynomial coefficients

};









#endif // CELESTIALBODY_H
