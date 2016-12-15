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
 *      160505    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */



#ifndef TAYLORSERIESINTEGRATION_H
#define TAYLORSERIESINTEGRATION_H


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

#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>
#include <Tudat/Mathematics/BasicMathematics/coordinateConversions.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h>

// Own classes

#include <thesisProject/celestialBody.h>            // Final version

#include <thesisProject/MarsAscentVehicle.h>    // Final version

#include <thesisProject/stateAndTime.h>             // Final version

#include <thesisProject/Auxiliary.h>                // Original test file

#include "thesisProject/projectLibraries/basicRecurrenceRelations.h"               // Original test file

#include "thesisProject/projectLibraries/allRecurrenceRelations.h"          // Original test file

#include <thesisProject/StepSize.h>             // Original test file

#include <thesisProject/projectLibraries/otherRequiredFunctions.h>      // deg2rad, rad2deg and B-P transformations



/// Main function ///
/// \brief performTaylorSeriesIntegrationStep is the primary integration function that is used to perform one integration step of TSI
/// \param planet
/// \param MAV
/// \param currentStateAndTime
/// \param stepSize
/// \param maxOrder
/// \return updatedStateAndTime This is a vector output with 8 elements: the updated three position values, the updated three velocity values, the updated mass value and the updated time.
///
Eigen::VectorXd performTaylorSeriesIntegrationStep(const celestialBody &planet_, const MarsAscentVehicle &MAV_, const StateAndTime &currentStateAndTime_, StepSize &stepSize, const double maxOrder_ = 20,
                                                   const double FlightPathAngle_ = 1000,
                                                   const double HeadingAngle_ = 1000);

// Ouput is the updated state and time

/* /// Some required functions /// (Moved to separate header and source file)


/// deg2rad ///
/// \brief deg2rad is a function to convert degrees to radians
/// \param deg
/// \return
///

const double deg2rad(const double deg);
//{

//    const double rad = deg*tudat::mathematical_constants::LONG_PI/180;

//    return rad;
//}

/// rad2deg ///
/// \brief rad2deg is a function to convert radians to degrees
/// \param rad
/// \return
///

const double rad2deg(const double rad);
//{

//    const double deg = rad*180/tudat::mathematical_constants::LONG_PI;

//    return deg;
//}


/// B-P frame transformations ///



//! Get transformation quaternion from the Body (B) to the Propulsion (P) frame.
Eigen::Quaterniond getBodyToPropulsionFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation );
//{
//    // Compute transformation quaternion.
//    // Note the sign change, because how angleAxisd is defined.
//    Eigen::AngleAxisd RotationAroundZaxis = Eigen::AngleAxisd(
//                -1.0 * thrustAzimuth, Eigen::Vector3d::UnitZ( ) );
//    Eigen::AngleAxisd RotationAroundYaxis = Eigen::AngleAxisd(
//                -1.0 * thrustElevation,
//                Eigen::Vector3d::UnitY( ) );
//    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
//                ( RotationAroundYaxis * RotationAroundZaxis ) );

//    // Return transformation quaternion.
//    return frameTransformationQuaternion;
//}


//! Get transformation matrix from the Body (B) to the Propulsion (P) frame.
Eigen::Matrix3d getBodyToPropulsionFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation );
//{
//    return getBodyToPropulsionFrameTransformationQuaternion(
//            thrustAzimuth, thrustElevation ).toRotationMatrix( );
//}

//! Get transformation matrix from the Propulsion (P) to the Body (B) frame.
Eigen::Matrix3d getPropulsionToBodyFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation );
//{
//    return getBodyToPropulsionFrameTransformationMatrix(
//            thrustAzimuth, thrustElevation ).transpose( );
//}

//! Get transformation quaternion from the Propulsion (P) to the Body (B) frame.
Eigen::Quaterniond getPropulsionToBodyFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation );
//{
//    return getBodyToPropulsionFrameTransformationQuaternion(
//            thrustAzimuth, thrustElevation ).inverse( );
//}

//*/

#endif // TAYLORSERIESINTEGRATION_H
