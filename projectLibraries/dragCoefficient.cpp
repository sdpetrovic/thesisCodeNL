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

#include "dragCoefficient.h"

namespace Drag
{

/// Drag coefficient function ///
/// \brief dragCoefficient  Computed the drag coefficient depending on the Mach number and the corresponding section of the C_D curve
/// \param machNumber   The current Mach number
/// \param dragCoefficientPolyCoefficients The drag coefficient polynomial curve fit coefficients
/// \param dragCoefficientMachRanges The drag coefficient section Mach ranges
/// \return
///

const double dragCoefficient(const double machNumber, const Eigen::MatrixXd dragCoefficientPolyCoefficients, const Eigen::MatrixXd dragCoefficientMachRanges){

    int section = 0;    // Defining the section and setting the default to 0

    for (int i = 0; i<6;i++){

        if (dragCoefficientMachRanges(i,0)<= machNumber && machNumber < dragCoefficientMachRanges(i,1)){

            section = i;
        }
    }


    const double currentDragCoefficient = dragCoefficientPolyCoefficients(section,1)*machNumber+dragCoefficientPolyCoefficients(section,0);

//std::cout<<"The CD section = "<<section<<std::endl;

return currentDragCoefficient;
} // end of dragCoefficient function

} // end of namespace Drag
