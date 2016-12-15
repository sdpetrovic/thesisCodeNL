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
 *      160824    S.D. Petrovic     Changed thrust angle dependency from time to altitude
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
#include <thesisProject/projectLibraries/otherRequiredFunctions.h>


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
//            Thrust_ = 5300.0;             //[N]           Taken from Trinidad et al. 2012
            Thrust_ = 5.3;             //[kN]           Taken from Trinidad et al. 2012
            thrustResetValue_ = Thrust_; // [kN]

            // Specific Impulse
            specificImpulse_ = 328.6;        //[s]       Taken from Trinidad et al. 2012

            // Reference Area
//            referenceArea_ = 0.091;          // [m^2]    Taken from Trinidad et al. 2012
            referenceArea_ = 9.1e-8;          // [km^2]    Taken from Trinidad et al. 2012

            // MAV GLOM
            MAVmass_ = 227.0;               // [kg]     Taken from Trinidad et al. 2012

            // Final altitude
            finalAltitude_ = 320.0;         // [km]     Taken from Trinidad et al. 2012


            // Drag Coefficient Polynomial Coefficients
                dragCoefficientPolyCoefficients_ = Eigen::MatrixXd::Zero(6,2);    // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve

                // Section 1: Mach range 0 to 0.5
                dragCoefficientPolyCoefficients_(0,0) = 0.2;

                // Section 2: Mach range 0.5 to 1
                dragCoefficientPolyCoefficients_(1,0) = -2.482534153247274e-16; //-2.483e-16;
                dragCoefficientPolyCoefficients_(1,1) = 0.4;

                // Section 3: Mach range 1 to 1.3
                dragCoefficientPolyCoefficients_(2,0) = -0.166666666666666; //-0.167;
                dragCoefficientPolyCoefficients_(2,1) = 0.566666666666666; //0.567;

                // Section 4: Mach range 1.3 to 2.5
                dragCoefficientPolyCoefficients_(3,0) = 0.754166666666666; //0.754;
                dragCoefficientPolyCoefficients_(3,1) = -0.1416666666666666; //-0.142;

                // Section 5: Mach range 2.5 to 4
                dragCoefficientPolyCoefficients_(4,0) = 0.566666666666667; //0.567;
                dragCoefficientPolyCoefficients_(4,1) = -0.06666666666667; //-0.0667;

                // Section 6: Mach range > 4
                dragCoefficientPolyCoefficients_(5,0) = 0.3;


            // Mach range per section for the mach-drag coefficient curve
                dragCoefficientMachRanges_ = Eigen::MatrixXd::Zero(6,2);

                // Section 1
                dragCoefficientMachRanges_(0,0) = 0.0;   // Lower bound
                dragCoefficientMachRanges_(0,1) = 0.5;   // Upper bound

                // Section 2
                dragCoefficientMachRanges_(1,0) = 0.5;   // Lower bound
                dragCoefficientMachRanges_(1,1) = 1.0;   // Upper bound

                // Section 3
                dragCoefficientMachRanges_(2,0) = 1.0;   // Lower bound
                dragCoefficientMachRanges_(2,1) = 1.3;   // Upper bound

                // Section 4
                dragCoefficientMachRanges_(3,0) = 1.3;   // Lower bound
                dragCoefficientMachRanges_(3,1) = 2.5;   // Upper bound

                // Section 5
                dragCoefficientMachRanges_(4,0) = 2.5;   // Lower bound
                dragCoefficientMachRanges_(4,1) = 4.0;   // Upper bound

                // Section 6
                dragCoefficientMachRanges_(5,0) = 4.0;   // Lower bound
                dragCoefficientMachRanges_(5,1) = 100.0; // Upper bound








        }

        else{

            // Setting all parameters to 0

            Thrust_ = 0.0;                                                    // T     engine nominal thrust
            specificImpulse_ = 0.0;                                           // Isp     engine nominal specific impulse
            referenceArea_ = 0.0;                                             // S   vehicle reference area
            MAVmass_ = 0.0;                                                   // m  vehicle GLOM
            finalAltitude_ = 0.0;                                             // h  final MAV altitude
            thrustResetValue_ = 0.0;                                          // Thrust reset value
            dragCoefficientPolyCoefficients_ = Eigen::MatrixXd::Zero(1,1);  // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
            dragCoefficientMachRanges_ = Eigen::MatrixXd::Zero(1,1);        // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient
//            thrustAzimuth_ = Eigen::MatrixXd::Zero(1,1);                    // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
//            thrustElevation_ = Eigen::MatrixXd::Zero(1,1);                  // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)




            std::cerr<<mav<<" is either not a valid string or has not yet been defined"<<std::endl;
        };


        /// Both thrust angles and the corresponding times will be optimised using the optimiser but can be set for validation purposes. Default = zeros ///

    // Thrust Azimuth-Gimbal Angles
        thrustAzimuth_ = Eigen::MatrixXd::Zero(6,3); // psiT   these are the thrust azimuth-gimbal angles in radians! as a function of altitude (including the altitude ranges)


        allTheSameAngleAzimuth_ = -0.299; /*0.7273;*/ // used in case they should all be the same

        // Section 1
        thrustAzimuth_(0,0) = -0.6;   // Lower bound altitude
        thrustAzimuth_(0,1) = -0.5;   // Upper bound altitude // 1.0

        thrustAzimuth_(0,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 2
        thrustAzimuth_(1,0) = thrustAzimuth_(0,1);   // Lower bound altitude
        thrustAzimuth_(1,1) = 53.66/2.0;   // Upper bound altitude // 5.0

        thrustAzimuth_(1,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 3
        thrustAzimuth_(2,0) = thrustAzimuth_(1,1);   // Lower bound altitude
        thrustAzimuth_(2,1) = 53.66;   // Upper bound altitude // 15.0

        thrustAzimuth_(2,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 4
        thrustAzimuth_(3,0) = thrustAzimuth_(2,1);   // Lower bound altitude
        thrustAzimuth_(3,1) = 60.0;   // Upper bound altitude // 35.0

        thrustAzimuth_(3,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 5
        thrustAzimuth_(4,0) = thrustAzimuth_(3,1);   // Lower bound altitude
        thrustAzimuth_(4,1) = 100.0;   // Upper bound altitude

        thrustAzimuth_(4,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 6
        thrustAzimuth_(5,0) = thrustAzimuth_(4,1);   // Lower bound altitude
        thrustAzimuth_(5,1) = finalAltitude_;   // Upper bound altitude

        thrustAzimuth_(5,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle


        // Thrust Elevation-Gimbal Angles
            thrustElevation_ = Eigen::MatrixXd::Zero(6,3); // epsilonT   these are the thrust elevation-gimbal angles in radians! as a function of altitude (including the altitude ranges)

            allTheSameAngleElevation_ = -0.184; /*-0.674;*/ // used in case they should all be the same

            // Section 1
            thrustElevation_(0,0) = -0.6;   // Lower bound altitude
            thrustElevation_(0,1) = -0.5;   // Upper bound altitude // 0.007

            thrustElevation_(0,2) = deg2rad(allTheSameAngleElevation_);   // Thrust elevation angle

            // Section 2
            thrustElevation_(1,0) = thrustElevation_(0,1);   // Lower bound altitude
            thrustElevation_(1,1) = 53.66/2.0;   // Upper bound altitude // 5.0

            thrustElevation_(1,2) = deg2rad(allTheSameAngleElevation_);   // Thrust elevation angle

            // Section 3
            thrustElevation_(2,0) = thrustElevation_(1,1);   // Lower bound altitude
            thrustElevation_(2,1) = 53.66;   // Upper bound altitude // 15.0

            thrustElevation_(2,2) = deg2rad(allTheSameAngleElevation_);   // Thrust elevation angle

            // Section 4
            thrustElevation_(3,0) = thrustElevation_(2,1);   // Lower bound altitude
            thrustElevation_(3,1) = 60;   // Upper bound altitude // 35.0

            thrustElevation_(3,2) = deg2rad(allTheSameAngleElevation_);   // Thrust elevation angle

            // Section 5
            thrustElevation_(4,0) = thrustElevation_(3,1);   // Lower bound altitude
            thrustElevation_(4,1) = 100.0;   // Upper bound altitude

            thrustElevation_(4,2) = deg2rad(allTheSameAngleElevation_);   // Thrust elevation angle

            // Section 6
            thrustElevation_(5,0) = thrustElevation_(4,1);   // Lower bound altitude
            thrustElevation_(5,1) = finalAltitude_;   // Upper bound altitude

            thrustElevation_(5,2) = deg2rad(allTheSameAngleElevation_);   // Thrust elevation angle


            /// Determine the thrust angle coefficients for the fit

            thrustAzimuthPolyCoefficients_ = Eigen::MatrixXd::Zero(thrustAzimuth_.rows(),2); // Set all polynomial coefficients to zero

            for (int l = 0;l<thrustAzimuth_.rows();l++){

                // Set the variables

                double x1,x2,y1,y2;


                if (l==0){
                    y1 = 0.0;
                }
                else {
                    y1 = thrustAzimuth_(l-1,2);
                }

                x1 = thrustAzimuth_(l,0);
                x2 = thrustAzimuth_(l,1);

                y2 = thrustAzimuth_(l,2);

                // Determine the coefficients (y = a*x+b)

                thrustAzimuthPolyCoefficients_(l,0) = (y2-y1)/(x2-x1); // a
                thrustAzimuthPolyCoefficients_(l,1) = (1/((1/x2)-(1/x1)))*((y2/x2)-(y1/x1)); // b

                /// Debug ///
//                std::cout<<"/// Debug Azimuth ///"<<std::endl;
//                std::cout<<"x1 = "<<x1<<std::endl;
//                std::cout<<"x2 = "<<x2<<std::endl;
//                std::cout<<"y1 = "<<y1<<std::endl;
//                std::cout<<"y2 = "<<y2<<std::endl;
//                std::cout<<"y2-y1 = "<<y2-y1<<std::endl;
//                std::cout<<"thrustAzimuthPolyCoefficients = "<<thrustAzimuthPolyCoefficients_<<std::endl;
//                std::cout<<"/// Debug Azimuth ///"<<std::endl;
                /// Debug ///

            } // End of set fit

            thrustElevationPolyCoefficients_ = Eigen::MatrixXd::Zero(thrustElevation_.rows(),2); // Set all polynomial coefficients to zero

            for (int l = 0;l<thrustElevation_.rows();l++){

                // Set the variables

                double x1,x2,y1,y2;


                if (l==0){
                    y1 = 0.0;
                }
                else {
                    y1 = thrustElevation_(l-1,2);
                }

                x1 = thrustElevation_(l,0);
                x2 = thrustElevation_(l,1);

                y2 = thrustElevation_(l,2);

                // Determine the coefficients (y = a*x+b)

                thrustElevationPolyCoefficients_(l,0) = (y2-y1)/(x2-x1); // a
                thrustElevationPolyCoefficients_(l,1) = (1/((1/x2)-(1/x1)))*((y2/x2)-(y1/x1)); // b

                /// Debug ///
//                std::cout<<"/// Debug Elevation ///"<<std::endl;
//                std::cout<<"x1 = "<<x1<<std::endl;
//                std::cout<<"x2 = "<<x2<<std::endl;
//                std::cout<<"y1 = "<<y1<<std::endl;
//                std::cout<<"y2 = "<<y2<<std::endl;
//                std::cout<<"y2-y1 = "<<y2-y1<<std::endl;
//                std::cout<<"thrustElevationPolyCoefficients = "<<thrustElevationPolyCoefficients_<<std::endl;
//                std::cout<<"/// Debug Elevation ///"<<std::endl;
                /// Debug ///

            } // End of set fit

                                                    } // End of constructor

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////// In order to add an extra MAV you simple add another elseif statement and provide the required information ////////////////
    ///////// PLEASE NOTE THAT IF NOT ALL INFORMATION IS PROVIDED IT COULD BE THAT A VALUE CLOSE TO 0 IS PROVIDED INSTEAD!//////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // Returning the different constant parameters

    const double Thrust() { return Thrust_; }                                                           // T     engine nominal thrust
    const double specificImpulse() { return specificImpulse_; }                                         // Isp     engine nominal specific impulse
    const double referenceArea() { return referenceArea_; }                                             // S   vehicle reference area
    const double MAVmass(){ return MAVmass_; }                                                          // m    vehicle GLOM

    // Returning the different polynomial coefficient parameter matrices

    const Eigen::MatrixXd dragCoefficientPolyCoefficients() { return dragCoefficientPolyCoefficients_; }            // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
    const Eigen::MatrixXd dragCoefficientMachRanges() { return dragCoefficientMachRanges_; }                        // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient


    // Returning the thrust angles as a function of time

    const Eigen::MatrixXd thrustAzimuth() { return thrustAzimuth_; }                                    // psiT   these are the thrust azimuth-gimbal angles as a function of time (including the time ranges)
    const Eigen::MatrixXd thrustElevation() { return thrustElevation_; }                                // epsilonT   these are the thrust elevation-gimbal angles as a function of time (including the time ranges)

    const Eigen::MatrixXd thrustAzimuthPolyCoefficient() { return thrustAzimuthPolyCoefficients_; }     // P_psiT   these are the polynomial coefficients for the fit for the azimuth angle curve
    const Eigen::MatrixXd thrustElevationPolyCoefficient() { return thrustElevationPolyCoefficients_; }     // P_epsilonT   these are the polynomial coefficients for the fit for the elevation angle curve



//// Set functions ////


    /// Thrust angles start ///
    void setThrustAzimuthSet(const Eigen::MatrixXd updatedThrustAzimuthSet)            // This function lets you provide the class with your own thrust azimuth angle set including altitude ranges
    {
        thrustAzimuth_ = updatedThrustAzimuthSet;

        // Set the coefficients for the fit
        thrustAzimuthPolyCoefficients_ = Eigen::MatrixXd::Zero(thrustAzimuth_.rows(),2); // Set all polynomial coefficient to zero

        for (int l = 0;l<thrustAzimuth_.rows();l++){

            // Set the variables

            double x1,x2,y1,y2;


            if (l==0){
                y1 = 0.0;
            }
            else {
                y1 = thrustAzimuth_(l-1,2);
            }

            x1 = thrustAzimuth_(l,0);
            x2 = thrustAzimuth_(l,1);

            y2 = thrustAzimuth_(l,2);

            // Determine the coefficients (y = a*x+b)

            thrustAzimuthPolyCoefficients_(l,0) = (y2-y1)/(x2-x1); // a
            thrustAzimuthPolyCoefficients_(l,1) = (1/((1/x2)-(1/x1)))*((y2/x2)-(y1/x1)); // b

            /// Debug ///
//            std::cout<<"/// Debug Azimuth ///"<<std::endl;
//            std::cout<<"x1 = "<<x1<<std::endl;
//            std::cout<<"x2 = "<<x2<<std::endl;
//            std::cout<<"y1 = "<<y1<<std::endl;
//            std::cout<<"y2 = "<<y2<<std::endl;
//            std::cout<<"thrustAzimuthPolyCoefficients = "<<thrustAzimuthPolyCoefficients_<<std::endl;
//            std::cout<<"/// Debug Azimuth ///"<<std::endl;
            /// Debug ///

        } // End of set fit
    }

    void setThrustElevationSet(const Eigen::MatrixXd updatedThrustElevationSet)            // This function lets you provide the class with your own thrust elevation angle set including altitude ranges
    {
        thrustElevation_ = updatedThrustElevationSet;

        // Set the coefficients for the fit
        thrustElevationPolyCoefficients_ = Eigen::MatrixXd::Zero(thrustElevation_.rows(),2); // Set all polynomial coefficients to zero

        for (int l = 0;l<thrustElevation_.rows();l++){

            // Set the variables

            double x1,x2,y1,y2;


            if (l==0){
                y1 = 0.0;
            }
            else {
                y1 = thrustElevation_(l-1,2);
            }

            x1 = thrustElevation_(l,0);
            x2 = thrustElevation_(l,1);

            y2 = thrustElevation_(l,2);

            // Determine the coefficients (y = a*x+b)

            thrustElevationPolyCoefficients_(l,0) = (y2-y1)/(x2-x1); // a
            thrustElevationPolyCoefficients_(l,1) = (1/((1/x2)-(1/x1)))*((y2/x2)-(y1/x1)); // b

            /// Debug ///
//            std::cout<<"/// Debug Elevation ///"<<std::endl;
//            std::cout<<"x1 = "<<x1<<std::endl;
//            std::cout<<"x2 = "<<x2<<std::endl;
//            std::cout<<"y1 = "<<y1<<std::endl;
//            std::cout<<"y2 = "<<y2<<std::endl;
//            std::cout<<"thrustElevationPolyCoefficients = "<<thrustElevationPolyCoefficients_<<std::endl;
//            std::cout<<"/// Debug Elevation ///"<<std::endl;
            /// Debug ///


        } // End of set fit
    }

    void setConstantThrustAzimuthAngle(const double updatedConstantThrustAzimuthAngle)      // This function lets you provide the class with your own constant thrust azimuth angle
    {
        allTheSameAngleAzimuth_ = updatedConstantThrustAzimuthAngle; // Update the constant angle and the values in the matrix

        // Section 1
        thrustAzimuth_(0,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 2
        thrustAzimuth_(1,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 3
        thrustAzimuth_(2,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 4
        thrustAzimuth_(3,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 5
        thrustAzimuth_(4,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

        // Section 6
        thrustAzimuth_(5,2) = deg2rad(allTheSameAngleAzimuth_);   // Thrust azimuth angle

    }


    void setConstantThrustElevationAngle(const double updatedConstantThrustElevationAngle)  // This function lets you provide the class with your own constan thrust elevation angle
    {
        allTheSameAngleElevation_ = updatedConstantThrustElevationAngle; // Update the constant angle and the values in the matrix

        // Section 1
        thrustElevation_(0,2) = deg2rad(allTheSameAngleElevation_);   // Thrust Elevation angle

        // Section 2
        thrustElevation_(1,2) = deg2rad(allTheSameAngleElevation_);   // Thrust Elevation angle

        // Section 3
        thrustElevation_(2,2) = deg2rad(allTheSameAngleElevation_);   // Thrust Elevation angle

        // Section 4
        thrustElevation_(3,2) = deg2rad(allTheSameAngleElevation_);   // Thrust Elevation angle

        // Section 5
        thrustElevation_(4,2) = deg2rad(allTheSameAngleElevation_);   // Thrust Elevation angle

        // Section 6
        thrustElevation_(5,2) = deg2rad(allTheSameAngleElevation_);   // Thrust Elevation angle

    }

    void setThrustAzimuthAngleValues(const Eigen::VectorXd updatedThrustAzimuthAngleVector)      // This function lets you provide the class with your own thrust azimuth angle values
    {

        // Section 1
        thrustAzimuth_(0,2) = deg2rad(updatedThrustAzimuthAngleVector(0));   // Thrust azimuth angle

        // Section 2
        thrustAzimuth_(1,2) = deg2rad(updatedThrustAzimuthAngleVector(1));   // Thrust azimuth angle

        // Section 3
        thrustAzimuth_(2,2) = deg2rad(updatedThrustAzimuthAngleVector(2));   // Thrust azimuth angle

        // Section 4
        thrustAzimuth_(3,2) = deg2rad(updatedThrustAzimuthAngleVector(3));   // Thrust azimuth angle

        // Section 5
        thrustAzimuth_(4,2) = deg2rad(updatedThrustAzimuthAngleVector(4));   // Thrust azimuth angle

        // Section 6
        thrustAzimuth_(5,2) = deg2rad(updatedThrustAzimuthAngleVector(5));   // Thrust azimuth angle

        // Set the coefficients for the fit
        thrustAzimuthPolyCoefficients_ = Eigen::MatrixXd::Zero(thrustAzimuth_.rows(),2); // Set all polynomial coefficient to zero

        for (int l = 0;l<thrustAzimuth_.rows();l++){

            // Set the variables

            double x1,x2,y1,y2;


            if (l==0){
                y1 = 0.0;
            }
            else {
                y1 = thrustAzimuth_(l-1,2);
            }

            x1 = thrustAzimuth_(l,0);
            x2 = thrustAzimuth_(l,1);

            y2 = thrustAzimuth_(l,2);

            // Determine the coefficients (y = a*x+b)

            thrustAzimuthPolyCoefficients_(l,0) = (y2-y1)/(x2-x1); // a
            thrustAzimuthPolyCoefficients_(l,1) = (1/((1/x2)-(1/x1)))*((y2/x2)-(y1/x1)); // b

            /// Debug ///
//            std::cout<<"/// Debug Azimuth ///"<<std::endl;
//            std::cout<<"x1 = "<<x1<<std::endl;
//            std::cout<<"x2 = "<<x2<<std::endl;
//            std::cout<<"y1 = "<<y1<<std::endl;
//            std::cout<<"y2 = "<<y2<<std::endl;
//            std::cout<<"thrustAzimuthPolyCoefficients = "<<thrustAzimuthPolyCoefficients_<<std::endl;
//            std::cout<<"/// Debug Azimuth ///"<<std::endl;
            /// Debug ///

        } // End of set fit
    }


    void setThrustElevationAngleValues(const Eigen::VectorXd updatedThrustElevationAngleVector)      // This function lets you provide the class with your own thrust Elevation angle values
    {

        // Section 1
        thrustElevation_(0,2) = deg2rad(updatedThrustElevationAngleVector(0));   // Thrust Elevation angle

        // Section 2
        thrustElevation_(1,2) = deg2rad(updatedThrustElevationAngleVector(1));   // Thrust Elevation angle

        // Section 3
        thrustElevation_(2,2) = deg2rad(updatedThrustElevationAngleVector(2));   // Thrust Elevation angle

        // Section 4
        thrustElevation_(3,2) = deg2rad(updatedThrustElevationAngleVector(3));   // Thrust Elevation angle

        // Section 5
        thrustElevation_(4,2) = deg2rad(updatedThrustElevationAngleVector(4));   // Thrust Elevation angle

        // Section 6
        thrustElevation_(5,2) = deg2rad(updatedThrustElevationAngleVector(5));   // Thrust Elevation angle

        // Set the coefficients for the fit
        thrustElevationPolyCoefficients_ = Eigen::MatrixXd::Zero(thrustElevation_.rows(),2); // Set all polynomial coefficients to zero

        for (int l = 0;l<thrustElevation_.rows();l++){

            // Set the variables

            double x1,x2,y1,y2;


            if (l==0){
                y1 = 0.0;
            }
            else {
                y1 = thrustElevation_(l-1,2);
            }

            x1 = thrustElevation_(l,0);
            x2 = thrustElevation_(l,1);

            y2 = thrustElevation_(l,2);

            // Determine the coefficients (y = a*x+b)

            thrustElevationPolyCoefficients_(l,0) = (y2-y1)/(x2-x1); // a
            thrustElevationPolyCoefficients_(l,1) = (1/((1/x2)-(1/x1)))*((y2/x2)-(y1/x1)); // b

            /// Debug ///
//            std::cout<<"/// Debug Elevation ///"<<std::endl;
//            std::cout<<"x1 = "<<x1<<std::endl;
//            std::cout<<"x2 = "<<x2<<std::endl;
//            std::cout<<"y1 = "<<y1<<std::endl;
//            std::cout<<"y2 = "<<y2<<std::endl;
//            std::cout<<"thrustElevationPolyCoefficients = "<<thrustElevationPolyCoefficients_<<std::endl;
//            std::cout<<"/// Debug Elevation ///"<<std::endl;
            /// Debug ///


        } // End of set fit
    }



    /// Thrust angles end set ///

    void setReferenceArea(const double updatedReferenceArea)                // This function can be used to change the reference area of the MAV. If set to 0 the drag can be neglected in all equations
    {
        referenceArea_ = updatedReferenceArea;
    }

    void setThrust(const double updatedThrust)                      // This function can be used to change the thrust value. If set to 0 the thrust can be neglected in all equations
    {
        Thrust_ = updatedThrust;
    }

    void setSpecificImpulse(const double updatedSpecificImpulse)        // This function can be used to change the specific impulse valuse.
    {
        specificImpulse_ = updatedSpecificImpulse;
    }

    void setMAVmass(const double updatedMass)               // This function can be used to change the MAV GLOM mass
    {
        MAVmass_ = updatedMass;
    }

    void setUpdatedFinalAltitude(const double updatedFinalAltitude)     // This function can be used to change the final altitude for the thrust elevation and azimuth ranges
    {
        int thrustAzimuthRow = thrustAzimuth_.rows()-1;
        int thrustElevationRow = thrustElevation_.rows()-1;

        thrustAzimuth_(thrustAzimuthRow,1) = updatedFinalAltitude;
        thrustElevation_(thrustElevationRow,1) = updatedFinalAltitude;

    }

    void setThrustResetValue(const double thrustResetValue){    // This function is used to make sure that then the thrust is reset, it is reset to the proper value!
        thrustResetValue_ = thrustResetValue;
    }

    //// Reset functions ////

    void resetThrust()
    {
        Thrust_ = thrustResetValue_ ; // kN                                   // This function is used to reset the thrust (theoretically)
    }

    void resetReferenceArea()
    {
        referenceArea_ = 9.1e-8; // [km^2]                                   // This function is used to reset the reference area, and thus drag (theoretically)
    }



protected:

private:

    // Creating the different constant parameters

    double Thrust_;                                                 // T     engine nominal thrust
    double specificImpulse_;                                        // Isp     engine nominal specific impulse
    double referenceArea_;                                          // S   vehicle reference area
    double MAVmass_;                                                // m    vehicle GLOM
    double finalAltitude_;                                          // h    final altitude of the MAV (pre-set)
    double thrustResetValue_;                                       // Thrust reset value
    double allTheSameAngleAzimuth_;                                 // psi  Set a constant thrust angle
    double allTheSameAngleElevation_;                               // epsilon    Set a constant thrust angle

    // Creating the different polynomial coefficient parameter matrices

    Eigen::MatrixXd dragCoefficientPolyCoefficients_;               // P_CDn     these are the polynomial coefficients for the fit for the drag coefficient curve
    Eigen::MatrixXd dragCoefficientMachRanges_;                     // dragCoefficientMachRanges      these are the Mach ranges corresponding to the polynomial coefficients for the drag coefficient

    // Creating the thrust angles as a function of time

    Eigen::MatrixXd thrustAzimuth_;                                 // psiT   these are the thrust azimuth-gimbal angles as a function of altitude (including the altitude ranges)
    Eigen::MatrixXd thrustElevation_;                               // epsilonT   these are the thrust elevation-gimbal angles as a function of altitude (including the altitude ranges)
    Eigen::MatrixXd thrustAzimuthPolyCoefficients_;                 // P_psiT   these are the polynomial coefficients for the fit for the psiT curve
    Eigen::MatrixXd thrustElevationPolyCoefficients_;               // P_epsilonT   these are the polynomial coefficients for the fit for the epsilonT curve



};









#endif // MARSASCENTVEHCILE_H
