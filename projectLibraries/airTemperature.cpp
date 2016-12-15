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

#include "airTemperature.h"

namespace air_temperature
{

const double airTemperature(const Eigen::MatrixXd temperaturePolyCoefficients, const Eigen::MatrixXd temperatureAltitudeRanges, double altitude){

    if (altitude <= -0.6){
        std::cerr<<"You have gone through the planet! The temperature will be set equal to the surface temperature"<<std::endl;
        altitude = -0.6;
    }

int section = 0;    // Define the section and setting the default to 0

for (int i = 0; i<5; i++){      // Determine in which section of the curve the MAV is

    if (temperatureAltitudeRanges(i,0)<=altitude && altitude < temperatureAltitudeRanges(i,1)){
        section = i;
    };
};

int order = 0; // Define the order and setting the default to 0

// Determine the order
if (section == 0){
    order = 1;

}
else if (section == 1){
    order = 2;
}
else if (section == 2){
    order = 6;
}
else if(section == 3){
    order = 8;
}
else if (section == 4){
    order = 0;
};


// Compute the current temperature

double currentTemperature = 0.0; // Set current temperature to 0;

for (int j = 0; j<order+1; j++){

    currentTemperature += temperaturePolyCoefficients(section,j)*pow(altitude,j);
};

//std::cout<<"The temperature section = "<<section<<std::endl;

return currentTemperature;

} // end of function airTemperature



} // end namespace air_temperature
