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
 *      160411    S.D. Petrovic          File created
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef MARSASCENTVEHICLE_H
#define MARSASCENTVEHICLE_H


#include <iostream>
#include <Eigen/Core>
#include <string>


// Author: Stacha Petrovic
// Date created: 11th April 2016


class MarsAscentVehicle

        /* This class will describe the different MAV characteristics
         * including thrust angle data.
         * These constants are:
         *
         * - T =                            [N]     engine nominal thrust
         * - Isp =                          [s]     engine nominal specific impulse
         * - S =                            [m^2]   vehicle reference area
         * - P_CDn                          [-]     these are the polynomial coefficients for the fit for the drag coefficient curve
         * - dragCoefficientMachRanges      [-]     these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient
         * - psiT                           [rad]   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
         * - epsilonT                       [rad]   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)
         *
         */





{
public:

    // Setting the celestial body with the default being the baseline liquid design (is the constructor)

    MarsAscentVehicle(const std::string mav = "LiteratureBaseline"){

        std::cout<<"You have chosen "<<mav<<" as your ascent vehicle"<<std::endl;

        if (mav == "LiteratureBaseline" || mav == "literatureBaseline" || mav == "literaturebaseline" || mav == "Baseline" || mav == "baseline")
        {

            /// Set different Baseline constants: ///


            // Thrust
            Thrust_ = 5300;             //[N]           Taken from Trinidad et al. 2012

            // Specific Impulse
            specificImpulse_ = 328.6;        //[s]       Taken from Trinidad et al. 2012

            // Reference Area
            referenceArea_ = 0.091;          // [m^2]    Taken from Trinidad et al. 2012

            // Drag Coefficient Polynomial Coefficients
                dragCoefficientPolyCoefficients_ = Eigen::MatrixXd::Zero(6,2);    // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve

                // Section 1: Mach range 0 to 0.5
                dragCoefficientPolyCoefficients_(0,0) = 0.2;

                // Section 2: Mach range 0.5 to 1
                dragCoefficientPolyCoefficients_(1,0) = -2.483e-16;
                dragCoefficientPolyCoefficients_(1,1) = 0.4;

                // Section 3: Mach range 1 to 1.3
                dragCoefficientPolyCoefficients_(2,0) = -0.167;
                dragCoefficientPolyCoefficients_(2,1) = 0.567;

                // Section 4: Mach range 1.3 to 2.5
                dragCoefficientPolyCoefficients_(3,0) = 0.754;
                dragCoefficientPolyCoefficients_(3,1) = -0.142;

                // Section 5: Mach range 2.5 to 4
                dragCoefficientPolyCoefficients_(4,0) = 0.567;
                dragCoefficientPolyCoefficients_(4,1) = -0.0667;

                // Section 6: Mach range > 4
                dragCoefficientPolyCoefficients_(5,0) = 0.3;


            // Mach range per section for the mach-drag coefficient curve
                dragCoefficientMachranges_ = Eigen::MatrixXd::Zero(6,2);

                // Section 1
                dragCoefficientMachranges_(0,0) = 0.0;   // Lower bound
                dragCoefficientMachranges_(0,1) = 0.5;   // Upper bound

                // Section 2
                dragCoefficientMachranges_(1,0) = 0.5;   // Lower bound
                dragCoefficientMachranges_(1,1) = 1.0;   // Upper bound

                // Section 3
                dragCoefficientMachranges_(2,0) = 1.0;   // Lower bound
                dragCoefficientMachranges_(2,1) = 1.3;   // Upper bound

                // Section 4
                dragCoefficientMachranges_(3,0) = 1.3;   // Lower bound
                dragCoefficientMachranges_(3,1) = 2.5;   // Upper bound

                // Section 5
                dragCoefficientMachranges_(4,0) = 2.5;   // Lower bound
                dragCoefficientMachranges_(4,1) = 4.0;   // Upper bound

                // Section 6
                dragCoefficientMachranges_(5,0) = 4.0;   // Lower bound
                dragCoefficientMachranges_(5,1) = 100.0; // Upper bound








        }

        else{

            // Setting all paramters to 0

            Thrust_ = 0;                                                    // T     engine nominal thrust
            specificImpulse_ = 0;                                           // Isp     engine nominal specific impulse
            referenceArea_ = 0;                                             // S   vehicle reference area
            dragCoefficientPolyCoefficients_ = Eigen::MatrixXd::Zero(1,1);  // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
            dragCoefficientMachranges_ = Eigen::MatrixXd::Zero(1,1);        // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient
//            thrustAzimuth_ = Eigen::MatrixXd::Zero(1,1);                    // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
//            thrustElevation_ = Eigen::MatrixXd::Zero(1,1);                  // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)




            std::cerr<<mav<<" is either not a valid string or has not yet been defined"<<std::endl;
        };


        /// Both thrust angles and the corresponding times will be optimised using the optimiser but can be set for validation purposes. Default = zeros ///

    // Thrust Azimuth-Gimbal Angles
        thrustAzimuth_ = Eigen::MatrixXd::Zero(6,3); // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)

        // Section 1
        thrustAzimuth_(0,0) = 0.0;   // Lower bound time
        thrustAzimuth_(0,1) = 0.0;   // Upper bound time

        thrustAzimuth_(0,2) = 0.0;   // Thrust azimuth angle

        // Section 2
        thrustAzimuth_(1,0) = 0.0;   // Lower bound time
        thrustAzimuth_(1,1) = 0.0;   // Upper bound time

        thrustAzimuth_(1,2) = 0.0;   // Thrust azimuth angle

        // Section 3
        thrustAzimuth_(2,0) = 0.0;   // Lower bound time
        thrustAzimuth_(2,1) = 0.0;   // Upper bound time

        thrustAzimuth_(2,2) = 0.0;   // Thrust azimuth angle

        // Section 4
        thrustAzimuth_(3,0) = 0.0;   // Lower bound time
        thrustAzimuth_(3,1) = 0.0;   // Upper bound time

        thrustAzimuth_(3,2) = 0.0;   // Thrust azimuth angle

        // Section 5
        thrustAzimuth_(4,0) = 0.0;   // Lower bound time
        thrustAzimuth_(4,1) = 0.0;   // Upper bound time

        thrustAzimuth_(4,2) = 0.0;   // Thrust azimuth angle

        // Section 6
        thrustAzimuth_(5,0) = 0.0;   // Lower bound time
        thrustAzimuth_(5,1) = 0.0;   // Upper bound time

        thrustAzimuth_(5,2) = 0.0;   // Thrust azimuth angle


        // Thrust Elevation-Gimbal Angles
            thrustElevation_ = Eigen::MatrixXd::Zero(6,3); // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)

            // Section 1
            thrustElevation_(0,0) = 0.0;   // Lower bound time
            thrustElevation_(0,1) = 0.0;   // Upper bound time

            thrustElevation_(0,2) = 0.0;   // Thrust elevation angle

            // Section 2
            thrustElevation_(1,0) = 0.0;   // Lower bound time
            thrustElevation_(1,1) = 0.0;   // Upper bound time

            thrustElevation_(1,2) = 0.0;   // Thrust elevation angle

            // Section 3
            thrustElevation_(2,0) = 0.0;   // Lower bound time
            thrustElevation_(2,1) = 0.0;   // Upper bound time

            thrustElevation_(2,2) = 0.0;   // Thrust elevation angle

            // Section 4
            thrustElevation_(3,0) = 0.0;   // Lower bound time
            thrustElevation_(3,1) = 0.0;   // Upper bound time

            thrustElevation_(3,2) = 0.0;   // Thrust elevation angle

            // Section 5
            thrustElevation_(4,0) = 0.0;   // Lower bound time
            thrustElevation_(4,1) = 0.0;   // Upper bound time

            thrustElevation_(4,2) = 0.0;   // Thrust elevation angle

            // Section 6
            thrustElevation_(5,0) = 0.0;   // Lower bound time
            thrustElevation_(5,1) = 0.0;   // Upper bound time

            thrustElevation_(5,2) = 0.0;   // Thrust elevation angle

                                                    } // End of constructor

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////// In order to add an extra MAV you simple add another elseif statement and provide the required information ////////////////
    ///////// PLEASE NOTE THAT IF NOT ALL INFORMATION IS PROVIDED IT COULD BE THAT A VALUE CLOSE TO 0 IS PROVIDED INSTEAD!//////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // Returning the different constant parameters

    const double Thrust() { return Thrust_; }                                                           // T     engine nominal thrust
    const double specificImpulse() { return specificImpulse_; }                                         // Isp     engine nominal specific impulse
    const double referenceArea() { return referenceArea_; }                                             // S   vehicle reference area

    // Returning the different polynomial coefficient parameter matrices

    const Eigen::MatrixXd dragCoefficientPolyCoefficients() { return dragCoefficientPolyCoefficients_; }            // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
    const Eigen::MatrixXd dragCoefficientMachranges() { return dragCoefficientMachranges_; }                        // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient


    // Returning the thrust angles as a function of time

    const Eigen::MatrixXd thrustAzimuth() { return thrustAzimuth_; }                                    // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
    const Eigen::MatrixXd thrustElevation() { return thrustElevation_; }                                // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)




    //// Set functions ////

    void setThrustAzimuth(const Eigen::MatrixXd updatedThrustAzimuthSet)            // This functions lets you provide the class with
    {
        thrustAzimuth_ = updatedThrustAzimuthSet;
    }

    void setThrustElevation(const Eigen::MatrixXd updatedThrustElevationSet)            // This functions lets you provide the class with
    {
        thrustElevation_ = updatedThrustElevationSet;
    }



protected:

private:

    // Creating the different constant parameters

    double Thrust_;                                                 // T     engine nominal thrust
    double specificImpulse_;                                        // Isp     engine nominal specific impulse
    double referenceArea_;                                          // S   vehicle reference area

    // Creating the different polynomial coefficient parameter matrices

    Eigen::MatrixXd dragCoefficientPolyCoefficients_;               // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
    Eigen::MatrixXd dragCoefficientMachranges_;                     // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient

    // Creating the thrust angles as a function of time

    Eigen::MatrixXd thrustAzimuth_;                                 // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
    Eigen::MatrixXd thrustElevation_;                               // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)



};









#endif // MARSASCENTVEHCILE_H
