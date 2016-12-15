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
 *      160428    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */


#ifndef ALLRECURRENCERELATIONS_H
#define ALLRECURRENCERELATIONS_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>

#include "tudatApplications/thesisProject/projectLibraries/allRecurrenceRelations.h"  // Not sure if such a long path is needed, but just in case



/// Main function ///
/// \brief getTaylorCoefficients
/// \param planetCharacteristics
/// \param vehicleCharacteristics
/// \param thrustAccelerationsBframe
/// \param initialEquationsVector
/// \param initialDerivativesVector
/// \param initialFunctionsMatrix
/// \param currentTime
/// \return
///
Eigen::MatrixXd getTaylorCoefficients(const double adiabeticIndex_, const double specificGasConstant_, const double standardGravitationalParameter_, const double rotationalVelocity_, const double primeMeridianAngle_,
                      const double inertialFrameTime_, const double bodyReferenceRadius_, const Eigen::MatrixXd temperaturePolyCoefficients_, const Eigen::MatrixXd temperatureAltitudeRanges_,
                      const Eigen::VectorXd densityPolyCoefficients_, const double Thrust_, const double specificImpulse_,
                      const double referenceArea_, const Eigen::MatrixXd dragCoefficientPolyCoefficients_, const Eigen::MatrixXd dragCoefficientMachRanges_,
        const Eigen::VectorXd thrustAccelerationsBframe_,
        const Eigen::VectorXd initialEquationsVector_,
        const Eigen::VectorXd initialDerivativesVector_,
        const Eigen::MatrixXd initialFunctionsMatrix_,
        const double currentTime_,
        const int maxOrder_);


/// Declare the auxiliary Equation vectors/matrix ///

//Eigen::MatrixXd XMatrix;


/*
Eigen::VectorXd XVector1 = Eigen::VectorXd::Zero(maxOrder);     // X1
Eigen::VectorXd XVector2 = Eigen::VectorXd::Zero(maxOrder);     // X2
Eigen::VectorXd XVector3 = Eigen::VectorXd::Zero(maxOrder);     // X3
Eigen::VectorXd XVector4 = Eigen::VectorXd::Zero(maxOrder);     // X4
Eigen::VectorXd XVector5 = Eigen::VectorXd::Zero(maxOrder);     // X5
Eigen::VectorXd XVector6 = Eigen::VectorXd::Zero(maxOrder);     // X6
Eigen::VectorXd XVector7 = Eigen::VectorXd::Zero(maxOrder);     // X7
Eigen::VectorXd XVector8 = Eigen::VectorXd::Zero(maxOrder);     // X8
Eigen::VectorXd XVector9 = Eigen::VectorXd::Zero(maxOrder);     // X9
Eigen::VectorXd XVector10 = Eigen::VectorXd::Zero(maxOrder);     // X10

Eigen::VectorXd XVector11 = Eigen::VectorXd::Zero(maxOrder);     // X11
Eigen::VectorXd XVector12 = Eigen::VectorXd::Zero(maxOrder);     // X12
Eigen::VectorXd XVector13 = Eigen::VectorXd::Zero(maxOrder);     // X13
Eigen::VectorXd XVector14 = Eigen::VectorXd::Zero(maxOrder);     // X14
Eigen::VectorXd XVector15 = Eigen::VectorXd::Zero(maxOrder);     // X15
Eigen::VectorXd XVector16 = Eigen::VectorXd::Zero(maxOrder);     // X16
Eigen::VectorXd XVector17 = Eigen::VectorXd::Zero(maxOrder);     // X17
Eigen::VectorXd XVector18 = Eigen::VectorXd::Zero(maxOrder);     // X18
Eigen::VectorXd XVector19 = Eigen::VectorXd::Zero(maxOrder);     // X19
Eigen::VectorXd XVector20 = Eigen::VectorXd::Zero(maxOrder);     // X20

Eigen::VectorXd XVector21 = Eigen::VectorXd::Zero(maxOrder);     // X21
Eigen::VectorXd XVector22 = Eigen::VectorXd::Zero(maxOrder);     // X22
Eigen::VectorXd XVector23 = Eigen::VectorXd::Zero(maxOrder);     // X23
Eigen::VectorXd XVector24 = Eigen::VectorXd::Zero(maxOrder);     // X24
Eigen::VectorXd XVector25 = Eigen::VectorXd::Zero(maxOrder);     // X25
Eigen::VectorXd XVector26 = Eigen::VectorXd::Zero(maxOrder);     // X26
Eigen::VectorXd XVector27 = Eigen::VectorXd::Zero(maxOrder);     // X27
Eigen::VectorXd XVector28 = Eigen::VectorXd::Zero(maxOrder);     // X28
Eigen::VectorXd XVector29 = Eigen::VectorXd::Zero(maxOrder);     // X29
Eigen::VectorXd XVector30 = Eigen::VectorXd::Zero(maxOrder);     // X30

Eigen::VectorXd XVector31 = Eigen::VectorXd::Zero(maxOrder);     // X31
Eigen::VectorXd XVector32 = Eigen::VectorXd::Zero(maxOrder);     // X32
Eigen::VectorXd XVector33 = Eigen::VectorXd::Zero(maxOrder);     // X33
Eigen::VectorXd XVector34 = Eigen::VectorXd::Zero(maxOrder);     // X34
Eigen::VectorXd XVector35 = Eigen::VectorXd::Zero(maxOrder);     // X35
Eigen::VectorXd XVector36 = Eigen::VectorXd::Zero(maxOrder);     // X36
Eigen::VectorXd XVector37 = Eigen::VectorXd::Zero(maxOrder);     // X37
Eigen::VectorXd XVector38 = Eigen::VectorXd::Zero(maxOrder);     // X38
Eigen::VectorXd XVector39 = Eigen::VectorXd::Zero(maxOrder);     // X39
Eigen::VectorXd XVector40 = Eigen::VectorXd::Zero(maxOrder);     // X40

Eigen::VectorXd XVector41 = Eigen::VectorXd::Zero(maxOrder);     // X41
Eigen::VectorXd XVector42 = Eigen::VectorXd::Zero(maxOrder);     // X42
Eigen::VectorXd XVector43 = Eigen::VectorXd::Zero(maxOrder);     // X43
Eigen::VectorXd XVector44 = Eigen::VectorXd::Zero(maxOrder);     // X44
Eigen::VectorXd XVector45 = Eigen::VectorXd::Zero(maxOrder);     // X45
Eigen::VectorXd XVector46 = Eigen::VectorXd::Zero(maxOrder);     // X46
Eigen::VectorXd XVector47 = Eigen::VectorXd::Zero(maxOrder);     // X47
Eigen::VectorXd XVector48 = Eigen::VectorXd::Zero(maxOrder);     // X48

//*/


/// Declare the auxiliary Derivative vectors/matrix ///

//Eigen::MatrixXd UMatrix;

/*
Eigen::VectorXd UVector1 = Eigen::VectorXd::Zero(maxOrder);     // X1
Eigen::VectorXd UVector2 = Eigen::VectorXd::Zero(maxOrder);     // X2
Eigen::VectorXd UVector3 = Eigen::VectorXd::Zero(maxOrder);     // X3
Eigen::VectorXd UVector4 = Eigen::VectorXd::Zero(maxOrder);     // X4
Eigen::VectorXd UVector5 = Eigen::VectorXd::Zero(maxOrder);     // X5
Eigen::VectorXd UVector6 = Eigen::VectorXd::Zero(maxOrder);     // X6
Eigen::VectorXd UVector7 = Eigen::VectorXd::Zero(maxOrder);     // X7
Eigen::VectorXd UVector8 = Eigen::VectorXd::Zero(maxOrder);     // X8
Eigen::VectorXd UVector9 = Eigen::VectorXd::Zero(maxOrder);     // X9
Eigen::VectorXd UVector10 = Eigen::VectorXd::Zero(maxOrder);     // X10

Eigen::VectorXd UVector11 = Eigen::VectorXd::Zero(maxOrder);     // X11
Eigen::VectorXd UVector12 = Eigen::VectorXd::Zero(maxOrder);     // X12
Eigen::VectorXd UVector13 = Eigen::VectorXd::Zero(maxOrder);     // X13
Eigen::VectorXd UVector14 = Eigen::VectorXd::Zero(maxOrder);     // X14
Eigen::VectorXd UVector15 = Eigen::VectorXd::Zero(maxOrder);     // X15
Eigen::VectorXd UVector16 = Eigen::VectorXd::Zero(maxOrder);     // X16
Eigen::VectorXd UVector17 = Eigen::VectorXd::Zero(maxOrder);     // X17
Eigen::VectorXd UVector18 = Eigen::VectorXd::Zero(maxOrder);     // X18
Eigen::VectorXd UVector19 = Eigen::VectorXd::Zero(maxOrder);     // X19
Eigen::VectorXd UVector20 = Eigen::VectorXd::Zero(maxOrder);     // X20

Eigen::VectorXd UVector21 = Eigen::VectorXd::Zero(maxOrder);     // X21
Eigen::VectorXd UVector22 = Eigen::VectorXd::Zero(maxOrder);     // X22
Eigen::VectorXd UVector23 = Eigen::VectorXd::Zero(maxOrder);     // X23
Eigen::VectorXd UVector24 = Eigen::VectorXd::Zero(maxOrder);     // X24
Eigen::VectorXd UVector25 = Eigen::VectorXd::Zero(maxOrder);     // X25
Eigen::VectorXd UVector26 = Eigen::VectorXd::Zero(maxOrder);     // X26
Eigen::VectorXd UVector27 = Eigen::VectorXd::Zero(maxOrder);     // X27
Eigen::VectorXd UVector28 = Eigen::VectorXd::Zero(maxOrder);     // X28
Eigen::VectorXd UVector29 = Eigen::VectorXd::Zero(maxOrder);     // X29
Eigen::VectorXd UVector30 = Eigen::VectorXd::Zero(maxOrder);     // X30

Eigen::VectorXd UVector31 = Eigen::VectorXd::Zero(maxOrder);     // X31
Eigen::VectorXd UVector32 = Eigen::VectorXd::Zero(maxOrder);     // X32
Eigen::VectorXd UVector33 = Eigen::VectorXd::Zero(maxOrder);     // X33
Eigen::VectorXd UVector34 = Eigen::VectorXd::Zero(maxOrder);     // X34
Eigen::VectorXd UVector35 = Eigen::VectorXd::Zero(maxOrder);     // X35
Eigen::VectorXd UVector36 = Eigen::VectorXd::Zero(maxOrder);     // X36
Eigen::VectorXd UVector37 = Eigen::VectorXd::Zero(maxOrder);     // X37
Eigen::VectorXd UVector38 = Eigen::VectorXd::Zero(maxOrder);     // X38
Eigen::VectorXd UVector39 = Eigen::VectorXd::Zero(maxOrder);     // X39
Eigen::VectorXd UVector40 = Eigen::VectorXd::Zero(maxOrder);     // X40

Eigen::VectorXd UVector41 = Eigen::VectorXd::Zero(maxOrder);     // X41
Eigen::VectorXd UVector42 = Eigen::VectorXd::Zero(maxOrder);     // X42
Eigen::VectorXd UVector43 = Eigen::VectorXd::Zero(maxOrder);     // X43
Eigen::VectorXd UVector44 = Eigen::VectorXd::Zero(maxOrder);     // X44
Eigen::VectorXd UVector45 = Eigen::VectorXd::Zero(maxOrder);     // X45
Eigen::VectorXd UVector46 = Eigen::VectorXd::Zero(maxOrder);     // X46
Eigen::VectorXd UVector47 = Eigen::VectorXd::Zero(maxOrder);     // X47
Eigen::VectorXd UVector48 = Eigen::VectorXd::Zero(maxOrder);     // X48

//*/

/// Declare the auxiliary Function vectors ///

/*
Eigen::VectorXd WVector4_1;     // W4,1
Eigen::VectorXd WVector4_2;     // W4,2
Eigen::VectorXd WVector4_3;     // W4,3
Eigen::VectorXd WVector4_4;     // W4,4
Eigen::VectorXd WVector4_5;     // W4,5
Eigen::VectorXd WVector4_6;     // W4,6
Eigen::VectorXd WVector4_7;     // W4,7
Eigen::VectorXd WVector4_8;     // W4,8
Eigen::VectorXd WVector4_9;     // W4,9
Eigen::VectorXd WVector4_10;     // W4,10

Eigen::VectorXd WVector4_11;     // W4,11
Eigen::VectorXd WVector4_12;     // W4,12
Eigen::VectorXd WVector4_13;     // W4,13
Eigen::VectorXd WVector4_14;     // W4,14
Eigen::VectorXd WVector4_15;     // W4,15
Eigen::VectorXd WVector4_16;     // W4,16
Eigen::VectorXd WVector4_17;     // W4,17
Eigen::VectorXd WVector4_18;     // W4,18
Eigen::VectorXd WVector4_19;     // W4,19
Eigen::VectorXd WVector4_20;     // W4,20

Eigen::VectorXd WVector4_21;     // W4,21
Eigen::VectorXd WVector4_22;     // W4,22
Eigen::VectorXd WVector4_23;     // W4,23
Eigen::VectorXd WVector4_24;     // W4,24


Eigen::VectorXd WVector5_1;     // W5,1
Eigen::VectorXd WVector5_2;     // W5,2
Eigen::VectorXd WVector5_3;     // W5,3
Eigen::VectorXd WVector5_4;     // W5,4
Eigen::VectorXd WVector5_5;     // W5,5
Eigen::VectorXd WVector5_6;     // W5,6
Eigen::VectorXd WVector5_7;     // W5,7
Eigen::VectorXd WVector5_8;     // W5,8


Eigen::VectorXd WVector6_1;     // W6,1
Eigen::VectorXd WVector6_2;     // W6,2
Eigen::VectorXd WVector6_3;     // W6,3
Eigen::VectorXd WVector6_4;     // W6,4
Eigen::VectorXd WVector6_5;     // W6,5
Eigen::VectorXd WVector6_6;     // W6,6


Eigen::VectorXd WVector8_1;     // W8,1
Eigen::VectorXd WVector8_2;     // W8,2
Eigen::VectorXd WVector8_3;     // W8,3


Eigen::VectorXd WVector9;     // W9


Eigen::VectorXd WVector12_1;     // W12,1
Eigen::VectorXd WVector12_2;     // W12,2
Eigen::VectorXd WVector12_3;     // W12,3
Eigen::VectorXd WVector12_4;     // W12,4
Eigen::VectorXd WVector12_5;     // W12,5
Eigen::VectorXd WVector12_6;     // W12,6
Eigen::VectorXd WVector12_7;     // W12,7


Eigen::VectorXd WVector13_1;     // W13,1
Eigen::VectorXd WVector13_2;     // W13,2
Eigen::VectorXd WVector13_3;     // W13,3
Eigen::VectorXd WVector13_4;     // W13,4
Eigen::VectorXd WVector13_5;     // W13,5


Eigen::VectorXd WVector14_1;     // W14,1
Eigen::VectorXd WVector14_2;     // W14,2
Eigen::VectorXd WVector14_3;     // W14,3


Eigen::VectorXd WVector15_1;     // W15,1
Eigen::VectorXd WVector15_2;     // W15,2


Eigen::VectorXd WVector16;     // W16


Eigen::VectorXd WVector20;     // W20


Eigen::VectorXd WVector21_1;     // W21,1
Eigen::VectorXd WVector21_2;     // W21,2
Eigen::VectorXd WVector21_3;     // W21,3


Eigen::VectorXd WVector23_1;     // W23,1
Eigen::VectorXd WVector23_2;     // W23,2
Eigen::VectorXd WVector23_3;     // W23,3
Eigen::VectorXd WVector23_4;     // W23,4


Eigen::VectorXd WVector24_1;     // W24,1
Eigen::VectorXd WVector24_2;     // W24,2
Eigen::VectorXd WVector24_3;     // W24,3
Eigen::VectorXd WVector24_4;     // W24,4
Eigen::VectorXd WVector24_5;     // W24,5
Eigen::VectorXd WVector24_6;     // W24,6
Eigen::VectorXd WVector24_7;     // W24,7
Eigen::VectorXd WVector24_8;     // W24,8
Eigen::VectorXd WVector24_9;     // W24,9
Eigen::VectorXd WVector24_10;     // W24,10

Eigen::VectorXd WVector24_11;     // W24,11
Eigen::VectorXd WVector24_12;     // W24,12
Eigen::VectorXd WVector24_13;     // W24,13
Eigen::VectorXd WVector24_14;     // W24,14
Eigen::VectorXd WVector24_15;     // W24,15
Eigen::VectorXd WVector24_16;     // W24,16
Eigen::VectorXd WVector24_17;     // W24,17


Eigen::VectorXd WVector25_1;     // W25,1
Eigen::VectorXd WVector25_2;     // W25,2
Eigen::VectorXd WVector25_3;     // W25,3


Eigen::VectorXd WVector26_1;     // W26,1
Eigen::VectorXd WVector26_2;     // W26,2
Eigen::VectorXd WVector26_3;     // W26,3


Eigen::VectorXd WVector27_1;     // W27,1
Eigen::VectorXd WVector27_2;     // W27,2
Eigen::VectorXd WVector27_3;     // W27,3
Eigen::VectorXd WVector27_4;     // W27,4
Eigen::VectorXd WVector27_5;     // W27,5
Eigen::VectorXd WVector27_6;     // W27,6


Eigen::VectorXd WVector28;     // W28


Eigen::VectorXd WVector30_1;     // W30,1
Eigen::VectorXd WVector30_2;     // W30,2
Eigen::VectorXd WVector30_3;     // W30,3
Eigen::VectorXd WVector30_4;     // W30,4
Eigen::VectorXd WVector30_5;     // W30,5
Eigen::VectorXd WVector30_6;     // W30,6
Eigen::VectorXd WVector30_7;     // W30,7
Eigen::VectorXd WVector30_8;     // W30,8
Eigen::VectorXd WVector30_9;     // W30,9


Eigen::VectorXd WVector32_1;     // W32,1
Eigen::VectorXd WVector32_2;     // W32,2
Eigen::VectorXd WVector32_3;     // W32,3
Eigen::VectorXd WVector32_4;     // W32,4


Eigen::VectorXd WVector33;     // W33


Eigen::VectorXd WVector34_2;     // W34,2
Eigen::VectorXd WVector34_3;     // W34,3
Eigen::VectorXd WVector34_4;     // W34,4


Eigen::VectorXd WVector35_1;     // W35,1
Eigen::VectorXd WVector35_2;     // W35,2
Eigen::VectorXd WVector35_3;     // W35,3


Eigen::VectorXd WVector36;     // W36


Eigen::VectorXd WVector37_1;     // W37,1
Eigen::VectorXd WVector37_2;     // W37,2
Eigen::VectorXd WVector37_3;     // W37,3
Eigen::VectorXd WVector37_4;     // W37,4


Eigen::VectorXd WVector38_1;     // W38,1
Eigen::VectorXd WVector38_2;     // W38,2
Eigen::VectorXd WVector38_3;     // W38,3


Eigen::VectorXd WVector40_1;     // W40,1
Eigen::VectorXd WVector40_2;     // W40,2
Eigen::VectorXd WVector40_3;     // W40,3
Eigen::VectorXd WVector40_4;     // W40,4


Eigen::VectorXd WVector41_1;     // W41,1
Eigen::VectorXd WVector41_2;     // W41,2


Eigen::VectorXd WVector42_1;     // W42,1
Eigen::VectorXd WVector42_2;     // W42,2
Eigen::VectorXd WVector42_3;     // W42,3
Eigen::VectorXd WVector42_4;     // W42,4
Eigen::VectorXd WVector42_5;     // W42,5
Eigen::VectorXd WVector42_6;     // W42,6
Eigen::VectorXd WVector42_7;     // W42,7
Eigen::VectorXd WVector42_8;     // W42,8


Eigen::VectorXd WVector43_1;     // W43,1
Eigen::VectorXd WVector43_2;     // W43,2


Eigen::VectorXd WVector45_1;     // W45,1
Eigen::VectorXd WVector45_2;     // W45,2
Eigen::VectorXd WVector45_3;     // W45,3
Eigen::VectorXd WVector45_4;     // W45,4
Eigen::VectorXd WVector45_5;     // W45,5
Eigen::VectorXd WVector45_6;     // W45,6
Eigen::VectorXd WVector45_7;     // W45,7
Eigen::VectorXd WVector45_8;     // W45,8


Eigen::VectorXd WVector47_1;     // W47,1
Eigen::VectorXd WVector47_2;     // W47,2
Eigen::VectorXd WVector47_3;     // W47,3


Eigen::VectorXd WVector48_1;     // W48,1
Eigen::VectorXd WVector48_2;     // W48,2

//*/

/*
Eigen::VectorXd WVector4_1 = Eigen::VectorXd::Zero(maxOrder);     // W4,1
Eigen::VectorXd WVector4_2 = Eigen::VectorXd::Zero(maxOrder);     // W4,2
Eigen::VectorXd WVector4_3 = Eigen::VectorXd::Zero(maxOrder);     // W4,3
Eigen::VectorXd WVector4_4 = Eigen::VectorXd::Zero(maxOrder);     // W4,4
Eigen::VectorXd WVector4_5 = Eigen::VectorXd::Zero(maxOrder);     // W4,5
Eigen::VectorXd WVector4_6 = Eigen::VectorXd::Zero(maxOrder);     // W4,6
Eigen::VectorXd WVector4_7 = Eigen::VectorXd::Zero(maxOrder);     // W4,7
Eigen::VectorXd WVector4_8 = Eigen::VectorXd::Zero(maxOrder);     // W4,8
Eigen::VectorXd WVector4_9 = Eigen::VectorXd::Zero(maxOrder);     // W4,9
Eigen::VectorXd WVector4_10 = Eigen::VectorXd::Zero(maxOrder);     // W4,10

Eigen::VectorXd WVector4_11 = Eigen::VectorXd::Zero(maxOrder);     // W4,11
Eigen::VectorXd WVector4_12 = Eigen::VectorXd::Zero(maxOrder);     // W4,12
Eigen::VectorXd WVector4_13 = Eigen::VectorXd::Zero(maxOrder);     // W4,13
Eigen::VectorXd WVector4_14 = Eigen::VectorXd::Zero(maxOrder);     // W4,14
Eigen::VectorXd WVector4_15 = Eigen::VectorXd::Zero(maxOrder);     // W4,15
Eigen::VectorXd WVector4_16 = Eigen::VectorXd::Zero(maxOrder);     // W4,16
Eigen::VectorXd WVector4_17 = Eigen::VectorXd::Zero(maxOrder);     // W4,17
Eigen::VectorXd WVector4_18 = Eigen::VectorXd::Zero(maxOrder);     // W4,18
Eigen::VectorXd WVector4_19 = Eigen::VectorXd::Zero(maxOrder);     // W4,19
Eigen::VectorXd WVector4_20 = Eigen::VectorXd::Zero(maxOrder);     // W4,20

Eigen::VectorXd WVector4_21 = Eigen::VectorXd::Zero(maxOrder);     // W4,21
Eigen::VectorXd WVector4_22 = Eigen::VectorXd::Zero(maxOrder);     // W4,22
Eigen::VectorXd WVector4_23 = Eigen::VectorXd::Zero(maxOrder);     // W4,23
Eigen::VectorXd WVector4_24 = Eigen::VectorXd::Zero(maxOrder);     // W4,24


Eigen::VectorXd WVector5_1 = Eigen::VectorXd::Zero(maxOrder);     // W5,1
Eigen::VectorXd WVector5_2 = Eigen::VectorXd::Zero(maxOrder);     // W5,2
Eigen::VectorXd WVector5_3 = Eigen::VectorXd::Zero(maxOrder);     // W5,3
Eigen::VectorXd WVector5_4 = Eigen::VectorXd::Zero(maxOrder);     // W5,4
Eigen::VectorXd WVector5_5 = Eigen::VectorXd::Zero(maxOrder);     // W5,5
Eigen::VectorXd WVector5_6 = Eigen::VectorXd::Zero(maxOrder);     // W5,6
Eigen::VectorXd WVector5_7 = Eigen::VectorXd::Zero(maxOrder);     // W5,7
Eigen::VectorXd WVector5_8 = Eigen::VectorXd::Zero(maxOrder);     // W5,8


Eigen::VectorXd WVector6_1 = Eigen::VectorXd::Zero(maxOrder);     // W6,1
Eigen::VectorXd WVector6_2 = Eigen::VectorXd::Zero(maxOrder);     // W6,2
Eigen::VectorXd WVector6_3 = Eigen::VectorXd::Zero(maxOrder);     // W6,3
Eigen::VectorXd WVector6_4 = Eigen::VectorXd::Zero(maxOrder);     // W6,4
Eigen::VectorXd WVector6_5 = Eigen::VectorXd::Zero(maxOrder);     // W6,5
Eigen::VectorXd WVector6_6 = Eigen::VectorXd::Zero(maxOrder);     // W6,6


Eigen::VectorXd WVector8_1 = Eigen::VectorXd::Zero(maxOrder);     // W8,1
Eigen::VectorXd WVector8_2 = Eigen::VectorXd::Zero(maxOrder);     // W8,2
Eigen::VectorXd WVector8_3 = Eigen::VectorXd::Zero(maxOrder);     // W8,3


Eigen::VectorXd WVector9 = Eigen::VectorXd::Zero(maxOrder);     // W9


Eigen::VectorXd WVector12_1 = Eigen::VectorXd::Zero(maxOrder);     // W12,1
Eigen::VectorXd WVector12_2 = Eigen::VectorXd::Zero(maxOrder);     // W12,2
Eigen::VectorXd WVector12_3 = Eigen::VectorXd::Zero(maxOrder);     // W12,3
Eigen::VectorXd WVector12_4 = Eigen::VectorXd::Zero(maxOrder);     // W12,4
Eigen::VectorXd WVector12_5 = Eigen::VectorXd::Zero(maxOrder);     // W12,5
Eigen::VectorXd WVector12_6 = Eigen::VectorXd::Zero(maxOrder);     // W12,6
Eigen::VectorXd WVector12_7 = Eigen::VectorXd::Zero(maxOrder);     // W12,7


Eigen::VectorXd WVector13_1 = Eigen::VectorXd::Zero(maxOrder);     // W13,1
Eigen::VectorXd WVector13_2 = Eigen::VectorXd::Zero(maxOrder);     // W13,2
Eigen::VectorXd WVector13_3 = Eigen::VectorXd::Zero(maxOrder);     // W13,3
Eigen::VectorXd WVector13_4 = Eigen::VectorXd::Zero(maxOrder);     // W13,4
Eigen::VectorXd WVector13_5 = Eigen::VectorXd::Zero(maxOrder);     // W13,5


Eigen::VectorXd WVector14_1 = Eigen::VectorXd::Zero(maxOrder);     // W14,1
Eigen::VectorXd WVector14_2 = Eigen::VectorXd::Zero(maxOrder);     // W14,2
Eigen::VectorXd WVector14_3 = Eigen::VectorXd::Zero(maxOrder);     // W14,3


Eigen::VectorXd WVector15_1 = Eigen::VectorXd::Zero(maxOrder);     // W15,1
Eigen::VectorXd WVector15_2 = Eigen::VectorXd::Zero(maxOrder);     // W15,2


Eigen::VectorXd WVector16 = Eigen::VectorXd::Zero(maxOrder);     // W16


Eigen::VectorXd WVector20 = Eigen::VectorXd::Zero(maxOrder);     // W20


Eigen::VectorXd WVector21_1 = Eigen::VectorXd::Zero(maxOrder);     // W21,1
Eigen::VectorXd WVector21_2 = Eigen::VectorXd::Zero(maxOrder);     // W21,2
Eigen::VectorXd WVector21_3 = Eigen::VectorXd::Zero(maxOrder);     // W21,3


Eigen::VectorXd WVector23_1 = Eigen::VectorXd::Zero(maxOrder);     // W23,1
Eigen::VectorXd WVector23_2 = Eigen::VectorXd::Zero(maxOrder);     // W23,2
Eigen::VectorXd WVector23_3 = Eigen::VectorXd::Zero(maxOrder);     // W23,3
Eigen::VectorXd WVector23_4 = Eigen::VectorXd::Zero(maxOrder);     // W23,4


Eigen::VectorXd WVector24_1 = Eigen::VectorXd::Zero(maxOrder);     // W24,1
Eigen::VectorXd WVector24_2 = Eigen::VectorXd::Zero(maxOrder);     // W24,2
Eigen::VectorXd WVector24_3 = Eigen::VectorXd::Zero(maxOrder);     // W24,3
Eigen::VectorXd WVector24_4 = Eigen::VectorXd::Zero(maxOrder);     // W24,4
Eigen::VectorXd WVector24_5 = Eigen::VectorXd::Zero(maxOrder);     // W24,5
Eigen::VectorXd WVector24_6 = Eigen::VectorXd::Zero(maxOrder);     // W24,6
Eigen::VectorXd WVector24_7 = Eigen::VectorXd::Zero(maxOrder);     // W24,7
Eigen::VectorXd WVector24_8 = Eigen::VectorXd::Zero(maxOrder);     // W24,8
Eigen::VectorXd WVector24_9 = Eigen::VectorXd::Zero(maxOrder);     // W24,9
Eigen::VectorXd WVector24_10 = Eigen::VectorXd::Zero(maxOrder);     // W24,10

Eigen::VectorXd WVector24_11 = Eigen::VectorXd::Zero(maxOrder);     // W24,11
Eigen::VectorXd WVector24_12 = Eigen::VectorXd::Zero(maxOrder);     // W24,12
Eigen::VectorXd WVector24_13 = Eigen::VectorXd::Zero(maxOrder);     // W24,13
Eigen::VectorXd WVector24_14 = Eigen::VectorXd::Zero(maxOrder);     // W24,14
Eigen::VectorXd WVector24_15 = Eigen::VectorXd::Zero(maxOrder);     // W24,15
Eigen::VectorXd WVector24_16 = Eigen::VectorXd::Zero(maxOrder);     // W24,16
Eigen::VectorXd WVector24_17 = Eigen::VectorXd::Zero(maxOrder);     // W24,17


Eigen::VectorXd WVector25_1 = Eigen::VectorXd::Zero(maxOrder);     // W25,1
Eigen::VectorXd WVector25_2 = Eigen::VectorXd::Zero(maxOrder);     // W25,2
Eigen::VectorXd WVector25_3 = Eigen::VectorXd::Zero(maxOrder);     // W25,3


Eigen::VectorXd WVector26_1 = Eigen::VectorXd::Zero(maxOrder);     // W26,1
Eigen::VectorXd WVector26_2 = Eigen::VectorXd::Zero(maxOrder);     // W26,2
Eigen::VectorXd WVector26_3 = Eigen::VectorXd::Zero(maxOrder);     // W26,3


Eigen::VectorXd WVector27_1 = Eigen::VectorXd::Zero(maxOrder);     // W27,1
Eigen::VectorXd WVector27_2 = Eigen::VectorXd::Zero(maxOrder);     // W27,2
Eigen::VectorXd WVector27_3 = Eigen::VectorXd::Zero(maxOrder);     // W27,3
Eigen::VectorXd WVector27_4 = Eigen::VectorXd::Zero(maxOrder);     // W27,4
Eigen::VectorXd WVector27_5 = Eigen::VectorXd::Zero(maxOrder);     // W27,5
Eigen::VectorXd WVector27_6 = Eigen::VectorXd::Zero(maxOrder);     // W27,6


Eigen::VectorXd WVector28 = Eigen::VectorXd::Zero(maxOrder);     // W28


Eigen::VectorXd WVector30_1 = Eigen::VectorXd::Zero(maxOrder);     // W30,1
Eigen::VectorXd WVector30_2 = Eigen::VectorXd::Zero(maxOrder);     // W30,2
Eigen::VectorXd WVector30_3 = Eigen::VectorXd::Zero(maxOrder);     // W30,3
Eigen::VectorXd WVector30_4 = Eigen::VectorXd::Zero(maxOrder);     // W30,4
Eigen::VectorXd WVector30_5 = Eigen::VectorXd::Zero(maxOrder);     // W30,5
Eigen::VectorXd WVector30_6 = Eigen::VectorXd::Zero(maxOrder);     // W30,6
Eigen::VectorXd WVector30_7 = Eigen::VectorXd::Zero(maxOrder);     // W30,7
Eigen::VectorXd WVector30_8 = Eigen::VectorXd::Zero(maxOrder);     // W30,8
Eigen::VectorXd WVector30_9 = Eigen::VectorXd::Zero(maxOrder);     // W30,9


Eigen::VectorXd WVector32_1 = Eigen::VectorXd::Zero(maxOrder);     // W32,1
Eigen::VectorXd WVector32_2 = Eigen::VectorXd::Zero(maxOrder);     // W32,2
Eigen::VectorXd WVector32_3 = Eigen::VectorXd::Zero(maxOrder);     // W32,3
Eigen::VectorXd WVector32_4 = Eigen::VectorXd::Zero(maxOrder);     // W32,4


Eigen::VectorXd WVector33 = Eigen::VectorXd::Zero(maxOrder);     // W33


Eigen::VectorXd WVector34_2 = Eigen::VectorXd::Zero(maxOrder);     // W34,2
Eigen::VectorXd WVector34_3 = Eigen::VectorXd::Zero(maxOrder);     // W34,3
Eigen::VectorXd WVector34_4 = Eigen::VectorXd::Zero(maxOrder);     // W34,4


Eigen::VectorXd WVector35_1 = Eigen::VectorXd::Zero(maxOrder);     // W35,1
Eigen::VectorXd WVector35_2 = Eigen::VectorXd::Zero(maxOrder);     // W35,2
Eigen::VectorXd WVector35_3 = Eigen::VectorXd::Zero(maxOrder);     // W35,3


Eigen::VectorXd WVector36 = Eigen::VectorXd::Zero(maxOrder);     // W36


Eigen::VectorXd WVector37_1 = Eigen::VectorXd::Zero(maxOrder);     // W37,1
Eigen::VectorXd WVector37_2 = Eigen::VectorXd::Zero(maxOrder);     // W37,2
Eigen::VectorXd WVector37_3 = Eigen::VectorXd::Zero(maxOrder);     // W37,3
Eigen::VectorXd WVector37_4 = Eigen::VectorXd::Zero(maxOrder);     // W37,4


Eigen::VectorXd WVector38_1 = Eigen::VectorXd::Zero(maxOrder);     // W38,1
Eigen::VectorXd WVector38_2 = Eigen::VectorXd::Zero(maxOrder);     // W38,2
Eigen::VectorXd WVector38_3 = Eigen::VectorXd::Zero(maxOrder);     // W38,3


Eigen::VectorXd WVector40_1 = Eigen::VectorXd::Zero(maxOrder);     // W40,1
Eigen::VectorXd WVector40_2 = Eigen::VectorXd::Zero(maxOrder);     // W40,2
Eigen::VectorXd WVector40_3 = Eigen::VectorXd::Zero(maxOrder);     // W40,3
Eigen::VectorXd WVector40_4 = Eigen::VectorXd::Zero(maxOrder);     // W40,4


Eigen::VectorXd WVector41_1 = Eigen::VectorXd::Zero(maxOrder);     // W41,1
Eigen::VectorXd WVector41_2 = Eigen::VectorXd::Zero(maxOrder);     // W41,2


Eigen::VectorXd WVector42_1 = Eigen::VectorXd::Zero(maxOrder);     // W42,1
Eigen::VectorXd WVector42_2 = Eigen::VectorXd::Zero(maxOrder);     // W42,2
Eigen::VectorXd WVector42_3 = Eigen::VectorXd::Zero(maxOrder);     // W42,3
Eigen::VectorXd WVector42_4 = Eigen::VectorXd::Zero(maxOrder);     // W42,4
Eigen::VectorXd WVector42_5 = Eigen::VectorXd::Zero(maxOrder);     // W42,5
Eigen::VectorXd WVector42_6 = Eigen::VectorXd::Zero(maxOrder);     // W42,6
Eigen::VectorXd WVector42_7 = Eigen::VectorXd::Zero(maxOrder);     // W42,7
Eigen::VectorXd WVector42_8 = Eigen::VectorXd::Zero(maxOrder);     // W42,8


Eigen::VectorXd WVector43_1 = Eigen::VectorXd::Zero(maxOrder);     // W43,1
Eigen::VectorXd WVector43_2 = Eigen::VectorXd::Zero(maxOrder);     // W43,2


Eigen::VectorXd WVector45_1 = Eigen::VectorXd::Zero(maxOrder);     // W45,1
Eigen::VectorXd WVector45_2 = Eigen::VectorXd::Zero(maxOrder);     // W45,2
Eigen::VectorXd WVector45_3 = Eigen::VectorXd::Zero(maxOrder);     // W45,3
Eigen::VectorXd WVector45_4 = Eigen::VectorXd::Zero(maxOrder);     // W45,4
Eigen::VectorXd WVector45_5 = Eigen::VectorXd::Zero(maxOrder);     // W45,5
Eigen::VectorXd WVector45_6 = Eigen::VectorXd::Zero(maxOrder);     // W45,6
Eigen::VectorXd WVector45_7 = Eigen::VectorXd::Zero(maxOrder);     // W45,7
Eigen::VectorXd WVector45_8 = Eigen::VectorXd::Zero(maxOrder);     // W45,8


Eigen::VectorXd WVector47_1 = Eigen::VectorXd::Zero(maxOrder);     // W47,1
Eigen::VectorXd WVector47_2 = Eigen::VectorXd::Zero(maxOrder);     // W47,2
Eigen::VectorXd WVector47_3 = Eigen::VectorXd::Zero(maxOrder);     // W47,3


Eigen::VectorXd WVector48_1 = Eigen::VectorXd::Zero(maxOrder);     // W48,1
Eigen::VectorXd WVector48_2 = Eigen::VectorXd::Zero(maxOrder);     // W48,2
//*/

/// Declare the return matrix for x1 till x7 ///

//Eigen::MatrixXd stateTaylorCoefficients;

/// Extra declarations ///

//int sectionCD;                  // The CD section required

//Eigen::VectorXd W9IntermediateVector;   // Intermediate vector created to be able to use the basic recurrence relations (should have done this from the start...)

//Eigen::VectorXd onesVector;     // Vector containing ones where needed

                                // Extra polynomial coefficient vectors for addition of vectors
//Eigen::VectorXd densityPolyCoefficient_1;
//Eigen::VectorXd temperaturePolyCoefficient_1_2;
//Eigen::VectorXd temperaturePolyCoefficient_1_3;
//Eigen::VectorXd temperaturePolyCoefficient_1_4;







#endif // ALLRECURRENCERELATIONS_H
