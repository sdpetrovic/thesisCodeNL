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
 *      160517    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef ASCENTDRAGFORCE_H
#define ASCENTDRAGFORCE_H

#include <iostream>
#include <Eigen/Core>
#include <cmath>

#include <thesisProject/projectLibraries/dragCoefficient.h>     // Drag coefficient function

namespace Drag
{

/// Drag force function ///
/// \brief ascentDragForce Function that computes the drag in kN
/// \param rotationalVelocity       The velocity of the MAV in the rotational reference frame
/// \param temp         The current temperature
/// \param adiabeticIndex   The adiabetic index gamma_a of the Martian atmosphere
/// \param specificGasConstant The specific gas constant R_a_star of the Martian atmosphere
/// \param dragCoefficientPolyCoefficients The polynomial coefficients for the drag coefficient curves
/// \param dragCoefficientMachRanges    The Mach ranges corresponding to the different sections of the drag coefficient curves
/// \param referenceArea    The reference area S [km^2]
/// \param currentAirDensity   The current air density
/// \return
///

const double ascentDragForce(const double rotationalVelocity, const double temp, const double adiabeticIndex, const double specificGasConstant,
                             const Eigen::MatrixXd dragCoefficientPolyCoefficients, const Eigen::MatrixXd dragCoefficientMachRanges, const double referenceArea, const double currentAirDensity);

} // end namespace Drag

#endif // ASCENTDRAGFORCE_H
