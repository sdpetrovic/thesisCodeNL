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
 *      160408    S.D. Petrovic          File created
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
// Date created: 8th April 2016


class celestialBody

        /* This class will describe the different celestial body characteristics
         * including atmospheric data.
         * These constants are:
         *
         * - gamma_a = 1.35                                         adiabetic index for Mars
         * - Rstar = Ra/Ma =~ 1918.42635            [m^2/(s^2*K)]   specific gas constant for Mars (Ra = 8.3144598 m^2*kg/(s^2*K*mol) and Ma = 0.004334 kg/mol)
         * - PTn                                                    these are the polynomial coefficients for the fit for the temperature curves
         * - Prho n                                                 these are the polynomial coefficients for the fit for the density curves
         * - mu_M = 4.283*10^13                     [m^3/s^2]       standard gravitational parameter for Mars
         * - marsRotationalVelocity = 7.088*10^-5   [rad/s]         rotational velocity of Mars
         * - OmegaP = 0                             [rad]           relative angle between the prime meridian and the x-axis (for now it is chosen such that
         *                                                          the inertial frame was set at time t = start time and going through the prime meridian)
         * - t0 = 0                                 [s]             time between the start time and the time that the inertial frame was set
         *
         */





{
public:

    celestialBody(){}

    // Creating the different constant parameters

    const double adiabeticIndex();                        // gamma_a      adiabetic index
    const double specificGasConstant();                   // Rstar        specific gas constant
    const double standardGravitationalParameter();        // mu_M         standard gravitational parameter
    const double marsRotationalVelocity();                // rotational velocity of Mars
    const double primeMeridianAngle();                    // OmegaP       relative angle between the prime meridian and the x-axis
    const double inertialFrameTime();                     // t0           time between the start time and the time that the inertial frame was set




    // Creating the different polynomial coefficient parameter matrices

    const Eigen::MatrixXd temperaturePolyCoefficients();   // PTn    Temperature polynomial coefficients

    const Eigen::VectorXd densityPolyCoefficients();       // Prho n Density polynomial coefficients


protected:

private:
};



    // Because these are constant values that should not be changed, all constants are defined outside of the class
    // Also, none of the calling 'functions' require any input

    // Set different constants:

    const double celestialBody::adiabeticIndex()
    {
        double adiabeticIndex_;

        adiabeticIndex_ = 1.35;

        return adiabeticIndex_;
    }


    const double celestialBody::specificGasConstant()
    {
        double specificGasConstant_;
        const double massOfMars = 0.004334;         // [kg/mol]
        const double gasConstant = 8.3144598;       // [m^2*kg/(s^2*K*mol)]


        specificGasConstant_ = gasConstant/massOfMars;  // [m^2/(s^2*K)]

        return specificGasConstant_;
    }


    const double celestialBody::standardGravitationalParameter()
    {
        double standardGravitationalParameter_;

        standardGravitationalParameter_ = 4.283e13;  // [m^3/s^2]

        return standardGravitationalParameter_;
    }


    const double celestialBody::marsRotationalVelocity()
    {
        double marsRotationalVelocity_;

        marsRotationalVelocity_ = 7.088e-5;  // [rad/s]

        return marsRotationalVelocity_;
    }


    const double celestialBody::primeMeridianAngle()
    {
        double primeMeridianAngle_;

        primeMeridianAngle_ = 0;  // [rad]

        return primeMeridianAngle_;
    }


    const double celestialBody::inertialFrameTime()
    {
        double inertialFrameTime_;

        inertialFrameTime_ = 0;  // [s]

        return inertialFrameTime_;
    }


    const Eigen::MatrixXd celestialBody::temperaturePolyCoefficients()
    {
        Eigen::MatrixXd temperaturePolyCoefficients_ = Eigen::MatrixXd::Zero(5, 9);  // PTn      Temperature polynomial coefficients


        // Section 1: -0.6 to 5.04 km MOLA
        temperaturePolyCoefficients_(0,0) = 194.165;
        temperaturePolyCoefficients_(0,1) = 3.415;

        // Section 2: 5.04 to 35.53 km MOLA
        temperaturePolyCoefficients_(1,0) = 222.051;
        temperaturePolyCoefficients_(1,1) = -2.130;
        temperaturePolyCoefficients_(1,2) = 0.006;

        // Section 3: 35.53 to 75.07 km MOLA
        temperaturePolyCoefficients_(2,0) = -1.167e4;
        temperaturePolyCoefficients_(2,1) = 1.407e3;
        temperaturePolyCoefficients_(2,2) = -68.294;
        temperaturePolyCoefficients_(2,3) = 1.733;
        temperaturePolyCoefficients_(2,4) = -0.0243;
        temperaturePolyCoefficients_(2,5) = 1.785e-4;
        temperaturePolyCoefficients_(2,6) = -5.388e-7;

        // Section 4: 75.07 to 170.05 km MOLA
        temperaturePolyCoefficients_(3,0) = 2.236e5;
        temperaturePolyCoefficients_(3,1) = -1.523e4;
        temperaturePolyCoefficients_(3,2) = 447.378;
        temperaturePolyCoefficients_(3,3) = -7.405;
        temperaturePolyCoefficients_(3,4) = 0.076;
        temperaturePolyCoefficients_(3,5) = -4.862e-4;
        temperaturePolyCoefficients_(3,6) = 1.931e-6;
        temperaturePolyCoefficients_(3,7) = -4.328e-9;
        temperaturePolyCoefficients_(3,8) = 4.1942e-12;

        // Section 5: 170.05 to 320 km MOLA
        temperaturePolyCoefficients_(4,0) = 136.5;

        return temperaturePolyCoefficients_;
    }


    const Eigen::VectorXd celestialBody::densityPolyCoefficients()
    {
        Eigen::VectorXd densityPolyCoefficients_ = Eigen::VectorXd::Zero(11);  // Prho n Density polynomial coefficients

        densityPolyCoefficients_(0) = -4.172;
        densityPolyCoefficients_(1) = -0.0962;
        densityPolyCoefficients_(2) = 1.414e-3;
        densityPolyCoefficients_(3) = -9.604e-5;
        densityPolyCoefficients_(4) = 2.273e-6;
        densityPolyCoefficients_(5) = -2.884e-8;
        densityPolyCoefficients_(6) = 2.146e-10;
        densityPolyCoefficients_(7) = -9.620e-13;
        densityPolyCoefficients_(8) = 2.559e-15;
        densityPolyCoefficients_(9) = -3.724e-18;
        densityPolyCoefficients_(10) = 2.287e-21;



        return densityPolyCoefficients_;
    }



//    const Eigen::MatrixXd temperaturePolyCoefficients = Eigen::MatrixXd::Zero(4, 9);   // PTn      Temperature polynomial coefficients

//    const Eigen::VectorXd densityPolyCoefficients = Eigen::VectorXd::Zero(11);         // Prho n Density polynomial coefficients








#endif // CELESTIALBODY_H
