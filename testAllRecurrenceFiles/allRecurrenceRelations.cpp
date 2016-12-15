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



#include "allRecurrenceRelations.h"
#include "basicRecurrenceRelations.h"



Eigen::MatrixXd getTaylorCoefficients(const double adiabeticIndex_, const double specificGasConstant_, const double standardGravitationalParameter_, const double rotationalVelocity_, const double primeMeridianAngle_,
                      const double inertialFrameTime_, const double bodyReferenceRadius_, const Eigen::MatrixXd temperaturePolyCoefficients_, const Eigen::MatrixXd temperatureAltitudeRanges_,
                      const Eigen::VectorXd densityPolyCoefficients_, const double Thrust_, const double specificImpulse_,
                      const double referenceArea_, const Eigen::MatrixXd dragCoefficientPolyCoefficients_, const Eigen::MatrixXd dragCoefficientMachRanges_,
        const Eigen::VectorXd thrustAccelerationsBframe_,
        const Eigen::VectorXd initialEquationsVector_,
        const Eigen::VectorXd initialDerivativesVector_,
        const Eigen::MatrixXd initialFunctionsMatrix_,
        const double currentTime_,
        const int maxOrder_){


    // Create new variables to use in this function and be deleted afterwards... I hope... Update: Nope...
    const double adiabeticIndex = adiabeticIndex_;
    const double specificGasConstant = specificGasConstant_;
    const double standardGravitationalParameter = standardGravitationalParameter_;
    const double rotationalVelocity = rotationalVelocity_;
//    const double primeMeridianAngle = primeMeridianAngle_;
//    const double inertialFrameTime = inertialFrameTime_;
//    const double bodyReferenceRadius = bodyReferenceRadius_;
    const Eigen::MatrixXd temperaturePolyCoefficients = temperaturePolyCoefficients_;
    const Eigen::MatrixXd temperatureAltitudeRanges = temperatureAltitudeRanges_;
    const Eigen::VectorXd densityPolyCoefficients = densityPolyCoefficients_;
//    const double Thrust = Thrust_;
//    const double specificImpulse = specificImpulse_;
    const double referenceArea = referenceArea_;
    const Eigen::MatrixXd dragCoefficientPolyCoefficients = dragCoefficientPolyCoefficients_;
    const Eigen::MatrixXd dragCoefficientMachRanges = dragCoefficientMachRanges_;
    const Eigen::VectorXd thrustAccelerationsBframe = thrustAccelerationsBframe_;
    const Eigen::VectorXd initialEquationsVector = initialEquationsVector_;
    const Eigen::VectorXd initialDerivativesVector = initialDerivativesVector_;
    const Eigen::MatrixXd initialFunctionsMatrix = initialFunctionsMatrix_;
//    const double currentTime = currentTime_;
    const int maxOrder = maxOrder_;







    /// Set initial values ///




    Eigen::MatrixXd XMatrix = Eigen::MatrixXd::Zero(initialEquationsVector.size(),maxOrder+1);        // Create the initial X and U matrices
    Eigen::MatrixXd UMatrix = Eigen::MatrixXd::Zero(initialDerivativesVector.size(),maxOrder);

    /// Debug ///
/*
    std::cout<<"initialEquationsVector = "<<initialEquationsVector<<std::endl;
    std::cout<<"initialDerivativesVector = "<<initialDerivativesVector<<std::endl;
    std::cout<<"initialEquationsVector.size() = "<<initialEquationsVector.size()<<std::endl;
    std::cout<<"XMatrix.col(0).size() = "<<XMatrix.col(0).size()<<std::endl;
    std::cout<<"initialEquationsVector(48) = "<<initialEquationsVector(48)<<std::endl;
    std::cout<<"initialDerivativesVector(48) = "<<initialDerivativesVector(48)<<std::endl;
//*/

    for (int i = 1; i < initialEquationsVector.size(); i++){                      // Fill the matrices

      XMatrix(i,0) = initialEquationsVector(i);
      XMatrix(i,1) = initialDerivativesVector(i);
      UMatrix(i,0) = initialDerivativesVector(i);




    };


    // Create the auxiliary functions vectors

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

    // And fill them

    WVector4_1(0) = initialFunctionsMatrix(4,1);     // W4,1
    WVector4_2(0) = initialFunctionsMatrix(4,2);     // W4,2
    WVector4_3(0) = initialFunctionsMatrix(4,3);     // W4,3
    WVector4_4(0) = initialFunctionsMatrix(4,4);     // W4,4
    WVector4_5(0) = initialFunctionsMatrix(4,5);     // W4,5
    WVector4_6(0) = initialFunctionsMatrix(4,6);     // W4,6
    WVector4_7(0) = initialFunctionsMatrix(4,7);     // W4,7
    WVector4_8(0) = initialFunctionsMatrix(4,8);     // W4,8
    WVector4_9(0) = initialFunctionsMatrix(4,9);     // W4,9
    WVector4_10(0) = initialFunctionsMatrix(4,10);     // W4,10

    WVector4_11(0) = initialFunctionsMatrix(4,11);     // W4,11
    WVector4_12(0) = initialFunctionsMatrix(4,12);     // W4,12
    WVector4_13(0) = initialFunctionsMatrix(4,13);     // W4,13
    WVector4_14(0) = initialFunctionsMatrix(4,14);     // W4,14
    WVector4_15(0) = initialFunctionsMatrix(4,15);     // W4,15
    WVector4_16(0) = initialFunctionsMatrix(4,16);     // W4,16
    WVector4_17(0) = initialFunctionsMatrix(4,17);     // W4,17
    WVector4_18(0) = initialFunctionsMatrix(4,18);     // W4,18
    WVector4_19(0) = initialFunctionsMatrix(4,19);     // W4,19
    WVector4_20(0) = initialFunctionsMatrix(4,20);     // W4,20

    WVector4_21(0) = initialFunctionsMatrix(4,21);     // W4,21
    WVector4_22(0) = initialFunctionsMatrix(4,22);     // W4,22
    WVector4_23(0) = initialFunctionsMatrix(4,23);     // W4,23
    WVector4_24(0) = initialFunctionsMatrix(4,24);     // W4,24


    WVector5_1(0) = initialFunctionsMatrix(5,1);     // W5,1
    WVector5_2(0) = initialFunctionsMatrix(5,2);     // W5,2
    WVector5_3(0) = initialFunctionsMatrix(5,3);     // W5,3
    WVector5_4(0) = initialFunctionsMatrix(5,4);     // W5,4
    WVector5_5(0) = initialFunctionsMatrix(5,5);     // W5,5
    WVector5_6(0) = initialFunctionsMatrix(5,6);     // W5,6
    WVector5_7(0) = initialFunctionsMatrix(5,7);     // W5,7
    WVector5_8(0) = initialFunctionsMatrix(5,8);     // W5,8


    WVector6_1(0) = initialFunctionsMatrix(6,1);     // W6,1
    WVector6_2(0) = initialFunctionsMatrix(6,2);     // W6,2
    WVector6_3(0) = initialFunctionsMatrix(6,3);     // W6,3
    WVector6_4(0) = initialFunctionsMatrix(6,4);     // W6,4
    WVector6_5(0) = initialFunctionsMatrix(6,5);     // W6,5
    WVector6_6(0) = initialFunctionsMatrix(6,6);     // W6,6


    WVector8_1(0) = initialFunctionsMatrix(8,1);     // W8,1
    WVector8_2(0) = initialFunctionsMatrix(8,2);     // W8,2
    WVector8_3(0) = initialFunctionsMatrix(8,3);     // W8,3


    WVector9(0) = initialFunctionsMatrix(9,1);     // W9


    WVector12_1(0) = initialFunctionsMatrix(12,1);     // W12,1
    WVector12_2(0) = initialFunctionsMatrix(12,2);     // W12,2
    WVector12_3(0) = initialFunctionsMatrix(12,3);     // W12,3
    WVector12_4(0) = initialFunctionsMatrix(12,4);     // W12,4
    WVector12_5(0) = initialFunctionsMatrix(12,5);     // W12,5
    WVector12_6(0) = initialFunctionsMatrix(12,6);     // W12,6
    WVector12_7(0) = initialFunctionsMatrix(12,7);     // W12,7


    WVector13_1(0) = initialFunctionsMatrix(13,1);     // W13,1
    WVector13_2(0) = initialFunctionsMatrix(13,2);     // W13,2
    WVector13_3(0) = initialFunctionsMatrix(13,3);     // W13,3
    WVector13_4(0) = initialFunctionsMatrix(13,4);     // W13,4
    WVector13_5(0) = initialFunctionsMatrix(13,5);     // W13,5


    WVector14_1(0) = initialFunctionsMatrix(14,1);     // W14,1
    WVector14_2(0) = initialFunctionsMatrix(14,2);     // W14,2
    WVector14_3(0) = initialFunctionsMatrix(14,3);     // W14,3


    WVector15_1(0) = initialFunctionsMatrix(15,1);     // W15,1
    WVector15_2(0) = initialFunctionsMatrix(15,2);     // W15,2


    WVector16(0) = initialFunctionsMatrix(16,1);     // W16


    WVector20(0) = initialFunctionsMatrix(20,1);     // W20


    WVector21_1(0) = initialFunctionsMatrix(21,1);     // W21,1
    WVector21_2(0) = initialFunctionsMatrix(21,2);     // W21,2
    WVector21_3(0) = initialFunctionsMatrix(21,3);     // W21,3


    WVector23_1(0) = initialFunctionsMatrix(23,1);     // W23,1
    WVector23_2(0) = initialFunctionsMatrix(23,2);     // W23,2
    WVector23_3(0) = initialFunctionsMatrix(23,3);     // W23,3
    WVector23_4(0) = initialFunctionsMatrix(23,4);     // W23,4


    WVector24_1(0) = initialFunctionsMatrix(24,1);     // W24,1
    WVector24_2(0) = initialFunctionsMatrix(24,2);     // W24,2
    WVector24_3(0) = initialFunctionsMatrix(24,3);     // W24,3
    WVector24_4(0) = initialFunctionsMatrix(24,4);     // W24,4
    WVector24_5(0) = initialFunctionsMatrix(24,5);     // W24,5
    WVector24_6(0) = initialFunctionsMatrix(24,6);     // W24,6
    WVector24_7(0) = initialFunctionsMatrix(24,7);     // W24,7
    WVector24_8(0) = initialFunctionsMatrix(24,8);     // W24,8
    WVector24_9(0) = initialFunctionsMatrix(24,9);     // W24,9
    WVector24_10(0) = initialFunctionsMatrix(24,10);     // W24,10

    WVector24_11(0) = initialFunctionsMatrix(24,11);     // W24,11
    WVector24_12(0) = initialFunctionsMatrix(24,12);     // W24,12
    WVector24_13(0) = initialFunctionsMatrix(24,13);     // W24,13
    WVector24_14(0) = initialFunctionsMatrix(24,14);     // W24,14
    WVector24_15(0) = initialFunctionsMatrix(24,15);     // W24,15
    WVector24_16(0) = initialFunctionsMatrix(24,16);     // W24,16
    WVector24_17(0) = initialFunctionsMatrix(24,17);     // W24,17


    WVector25_1(0) = initialFunctionsMatrix(25,1);     // W25,1
    WVector25_2(0) = initialFunctionsMatrix(25,2);     // W25,2
    WVector25_3(0) = initialFunctionsMatrix(25,3);     // W25,3


    WVector26_1(0) = initialFunctionsMatrix(26,1);     // W26,1
    WVector26_2(0) = initialFunctionsMatrix(26,2);     // W26,2
    WVector26_3(0) = initialFunctionsMatrix(26,3);     // W26,3


    WVector27_1(0) = initialFunctionsMatrix(27,1);     // W27,1
    WVector27_2(0) = initialFunctionsMatrix(27,2);     // W27,2
    WVector27_3(0) = initialFunctionsMatrix(27,3);     // W27,3
    WVector27_4(0) = initialFunctionsMatrix(27,4);     // W27,4
    WVector27_5(0) = initialFunctionsMatrix(27,5);     // W27,5
    WVector27_6(0) = initialFunctionsMatrix(27,6);     // W27,6


    WVector28(0) = initialFunctionsMatrix(28,1);     // W28


    WVector30_1(0) = initialFunctionsMatrix(30,1);     // W30,1
    WVector30_2(0) = initialFunctionsMatrix(30,2);     // W30,2
    WVector30_3(0) = initialFunctionsMatrix(30,3);     // W30,3
    WVector30_4(0) = initialFunctionsMatrix(30,4);     // W30,4
    WVector30_5(0) = initialFunctionsMatrix(30,5);     // W30,5
    WVector30_6(0) = initialFunctionsMatrix(30,6);     // W30,6
    WVector30_7(0) = initialFunctionsMatrix(30,7);     // W30,7
    WVector30_8(0) = initialFunctionsMatrix(30,8);     // W30,8
    WVector30_9(0) = initialFunctionsMatrix(30,9);     // W30,9


    WVector32_1(0) = initialFunctionsMatrix(32,1);     // W32,1
    WVector32_2(0) = initialFunctionsMatrix(32,2);     // W32,2
    WVector32_3(0) = initialFunctionsMatrix(32,3);     // W32,3
    WVector32_4(0) = initialFunctionsMatrix(32,4);     // W32,4


    WVector33(0) = initialFunctionsMatrix(33,1);     // W33


    WVector34_2(0) = initialFunctionsMatrix(34,1);     // W34,2
    WVector34_3(0) = initialFunctionsMatrix(34,2);     // W34,3
    WVector34_4(0) = initialFunctionsMatrix(34,3);     // W34,4


    WVector35_1(0) = initialFunctionsMatrix(35,1);     // W35,1
    WVector35_2(0) = initialFunctionsMatrix(35,2);     // W35,2
    WVector35_3(0) = initialFunctionsMatrix(35,3);     // W35,3


    WVector36(0) = initialFunctionsMatrix(36,1);     // W36


    WVector37_1(0) = initialFunctionsMatrix(37,1);     // W37,1
    WVector37_2(0) = initialFunctionsMatrix(37,2);     // W37,2
    WVector37_3(0) = initialFunctionsMatrix(37,3);     // W37,3
    WVector37_4(0) = initialFunctionsMatrix(37,4);     // W37,4


    WVector38_1(0) = initialFunctionsMatrix(38,1);     // W38,1
    WVector38_2(0) = initialFunctionsMatrix(38,2);     // W38,2
    WVector38_3(0) = initialFunctionsMatrix(38,3);     // W38,3


    WVector40_1(0) = initialFunctionsMatrix(40,1);     // W40,1
    WVector40_2(0) = initialFunctionsMatrix(40,2);     // W40,2
    WVector40_3(0) = initialFunctionsMatrix(40,3);     // W40,3
    WVector40_4(0) = initialFunctionsMatrix(40,4);     // W40,4


    WVector41_1(0) = initialFunctionsMatrix(41,1);     // W41,1
    WVector41_2(0) = initialFunctionsMatrix(41,2);     // W41,2


    WVector42_1(0) = initialFunctionsMatrix(42,1);     // W42,1
    WVector42_2(0) = initialFunctionsMatrix(42,2);     // W42,2
    WVector42_3(0) = initialFunctionsMatrix(42,3);     // W42,3
    WVector42_4(0) = initialFunctionsMatrix(42,4);     // W42,4
    WVector42_5(0) = initialFunctionsMatrix(42,5);     // W42,5
    WVector42_6(0) = initialFunctionsMatrix(42,6);     // W42,6
    WVector42_7(0) = initialFunctionsMatrix(42,7);     // W42,7
    WVector42_8(0) = initialFunctionsMatrix(42,8);     // W42,8



    WVector43_1(0) = initialFunctionsMatrix(43,1);     // W43,1
    WVector43_2(0) = initialFunctionsMatrix(43,2);     // W43,2


    WVector45_1(0) = initialFunctionsMatrix(45,1);     // W45,1
    WVector45_2(0) = initialFunctionsMatrix(45,2);     // W45,2
    WVector45_3(0) = initialFunctionsMatrix(45,3);     // W45,3
    WVector45_4(0) = initialFunctionsMatrix(45,4);     // W45,4
    WVector45_5(0) = initialFunctionsMatrix(45,5);     // W45,5
    WVector45_6(0) = initialFunctionsMatrix(45,6);     // W45,6
    WVector45_7(0) = initialFunctionsMatrix(45,7);     // W45,7
    WVector45_8(0) = initialFunctionsMatrix(45,8);     // W45,8


    WVector47_1(0) = initialFunctionsMatrix(47,1);     // W47,1
    WVector47_2(0) = initialFunctionsMatrix(47,2);     // W47,2
    WVector47_3(0) = initialFunctionsMatrix(47,3);     // W47,3


    WVector48_1(0) = initialFunctionsMatrix(48,1);     // W48,1
    WVector48_2(0) = initialFunctionsMatrix(48,2);     // W48,2

    ////// Declaring other Vectors and variables ///////


    Eigen::VectorXd W9IntermediateVector = Eigen::VectorXd::Zero(maxOrder);         // Intermediate vector created to be able to use the basic recurrence relations (should have done this from the start...)
    W9IntermediateVector(0) = XMatrix(9,0)*UMatrix(8,0);            // And setting the first value


    Eigen::VectorXd onesVector = Eigen::VectorXd::Zero(maxOrder);       // Vector used to represent the int 1 in cases where the entire vector has to be substracted from 1(for instance)

    Eigen::VectorXd densityPolyCoefficient_1 = Eigen::VectorXd::Zero(maxOrder);             // Vectors used to represent the polynomials when required for addition of vectors
    Eigen::VectorXd temperaturePolyCoefficient_1_2 = Eigen::VectorXd::Zero(maxOrder);
    Eigen::VectorXd temperaturePolyCoefficient_1_3 = Eigen::VectorXd::Zero(maxOrder);
    Eigen::VectorXd temperaturePolyCoefficient_1_4 = Eigen::VectorXd::Zero(maxOrder);

    int sectionCD = 0;      // Drag coefficient function section and setting the default to zero

    /// Debug ///

//    std::cout<<"So far so good 1"<<std::endl;


/// Recurrence computations ///

    // getMultiplicationRecurrenceRelation                  Multiplication
    // getDivisionRecurrenceRelation                        Division
    // getPowerRecurrenceRelation                           Power
    // getExponentialRecurrenceRelation                     Exponential
    // getCosineRecurrenceRelation                          Cosine
    // getSineRecurrenceRelation                            Sine

    int count = 1;  // Counter

    for (int k = 1; k< maxOrder;k++){

//    std::cout<<"So far so good 2"<<std::endl;

        // Set general required vectors

        for (int i = 0; i< k+1; i++){

            onesVector(i) = 1;
            densityPolyCoefficient_1(i) = densityPolyCoefficients(1);
            temperaturePolyCoefficient_1_2(i) = temperaturePolyCoefficients(1,0);
            temperaturePolyCoefficient_1_3(i) = temperaturePolyCoefficients(2,0);
            temperaturePolyCoefficient_1_4(i) = temperaturePolyCoefficients(3,0);


        }


//    std::cout<<"So far so good 3"<<std::endl;

        // 1
        UMatrix(1,k) = XMatrix(4,k);

        // 2
        UMatrix(2,k) = XMatrix(5,k);

        // 3
        UMatrix(3,k) = XMatrix(6,k);

        // 4

        WVector4_1(k) = getDivisionRecurrenceRelation(XMatrix.row(1),XMatrix.row(9),WVector4_1,k); //  XMatrix(1,k)/XMatrix(9,k);
        WVector4_2(k) = thrustAccelerationsBframe(0)-getDivisionRecurrenceRelation(XMatrix.row(27),XMatrix.row(7),WVector4_2,k);
        WVector4_3(k) = getCosineRecurrenceRelation((XMatrix.row(10)+XMatrix.row(11)),WVector4_8,k);
        WVector4_4(k) = getSineRecurrenceRelation(XMatrix.row(12),WVector4_6,k);
        WVector4_5(k) = getCosineRecurrenceRelation(XMatrix.row(13),WVector4_9,k);
        WVector4_6(k) = getCosineRecurrenceRelation(XMatrix.row(12),WVector4_4,k);
        WVector4_7(k) = getSineRecurrenceRelation(XMatrix.row(14),XMatrix.row(16),k);
        WVector4_8(k) = getSineRecurrenceRelation((XMatrix.row(10)+XMatrix.row(11)),WVector4_3,k);
        WVector4_9(k) = getSineRecurrenceRelation(XMatrix.row(13),WVector4_5,k);
        WVector4_10(k) = getMultiplicationRecurrenceRelation(WVector4_4,WVector4_5,k);
        WVector4_11(k) = getMultiplicationRecurrenceRelation(WVector4_6,WVector4_7,k);
        WVector4_12(k) = getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(16),k);
        WVector4_13(k) = getMultiplicationRecurrenceRelation(WVector4_4,WVector4_9,k);
        WVector4_14(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_5,k);
        WVector4_15(k) = getMultiplicationRecurrenceRelation(WVector4_6,XMatrix.row(16),k);
        WVector4_16(k) = getMultiplicationRecurrenceRelation(WVector4_9,WVector4_7,k);
        WVector4_17(k) = getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(16),k);
        WVector4_18(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_12,k);
        WVector4_19(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_13,k);
        WVector4_20(k) = getMultiplicationRecurrenceRelation(WVector4_10,WVector4_7,k);
        WVector4_21(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_16,k);
        WVector4_22(k) = getMultiplicationRecurrenceRelation(WVector4_3,(WVector4_11-WVector4_17),k);
        WVector4_23(k) = getMultiplicationRecurrenceRelation(WVector4_3,(-WVector4_20-WVector4_15),k);
        WVector4_24(k) = getMultiplicationRecurrenceRelation(WVector4_2,(WVector4_22-WVector4_18),k);

        UMatrix(4,k) = -standardGravitationalParameter*WVector4_1(k)+WVector4_24(k)+thrustAccelerationsBframe(1)*(WVector4_19(k)-WVector4_14(k))+thrustAccelerationsBframe(2)*(WVector4_23(k)-WVector4_21(k));

//        std::cout<<"So far so good 4"<<std::endl;

        // 5

        WVector5_1(k) = getDivisionRecurrenceRelation(XMatrix.row(2),XMatrix.row(9),WVector5_1,k);
        WVector5_2(k) = getMultiplicationRecurrenceRelation(WVector4_8,(WVector4_11-WVector4_17),k);
        WVector5_3(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_12,k);
        WVector5_4(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_13,k);
        WVector5_5(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_5,k);
        WVector5_6(k) = getMultiplicationRecurrenceRelation(WVector4_8,(-WVector4_20-WVector4_11),k);
        WVector5_7(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_16,k);
        WVector5_8(k) = getMultiplicationRecurrenceRelation(WVector4_2,(WVector5_2+WVector5_3),k);


        UMatrix(5,k) = -standardGravitationalParameter*WVector5_1(k)+WVector5_8(k)+thrustAccelerationsBframe(1)*(WVector5_4(k)+WVector5_5(k))+thrustAccelerationsBframe(2)*(WVector5_6(k)+WVector5_7(k));

        // 6

        WVector6_1(k) = getDivisionRecurrenceRelation(XMatrix.row(3),XMatrix.row(9),WVector6_1,k);
        WVector6_2(k) = getMultiplicationRecurrenceRelation(WVector4_5,WVector4_15,k);
        WVector6_3(k) = getMultiplicationRecurrenceRelation(WVector4_6,WVector4_9,k);
        WVector6_4(k) = getMultiplicationRecurrenceRelation(WVector4_5,WVector4_11,k);
        WVector6_5(k) = getMultiplicationRecurrenceRelation(WVector4_4,XMatrix.row(16),k);
        WVector6_6(k) = getMultiplicationRecurrenceRelation(WVector4_2,(WVector6_2-WVector4_11),k);

        UMatrix(6,k) = -standardGravitationalParameter*WVector6_1(k)+WVector6_6(k)-thrustAccelerationsBframe(1)*WVector6_3(k)+thrustAccelerationsBframe(2)*(WVector6_4(k)-WVector6_5(k));

        // 7

//        UMatrix(7,k) =-thrust/(tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION*specificImpulse);
        UMatrix(7,k) = 0;

        // 8

        WVector8_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(1),XMatrix.row(4),k);
        WVector8_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(2),XMatrix.row(5),k);
        WVector8_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(3),XMatrix.row(6),k);

        UMatrix(8,k) = 2*(WVector8_1(k)+WVector8_2(k)+WVector8_3(k));

        // 10

        UMatrix(10,k) = 0;

        // 11

        UMatrix(11,k) = XMatrix(45,k)-rotationalVelocity;


        // 12

        WVector12_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(20),XMatrix.row(6),k);
        WVector12_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(3),XMatrix.row(25),k);
        WVector12_3(k) = getDivisionRecurrenceRelation(XMatrix.row(3),XMatrix.row(20),WVector12_3,k);
        WVector12_4(k) = getMultiplicationRecurrenceRelation(WVector12_3,WVector12_3,k);
        WVector12_5(k) = getPowerRecurrenceRelation((onesVector-WVector12_4),WVector12_5,0.5,k);                                 // change to Eigen::VectorXd -- check
        WVector12_6(k) = getMultiplicationRecurrenceRelation(XMatrix.row(8),WVector12_5,k);
        WVector12_7(k) = getDivisionRecurrenceRelation((WVector12_1-WVector12_2),WVector12_6,WVector12_7,k);


        UMatrix(12,k) = WVector12_7(k);

        // 19

        UMatrix(19,k) = 2*(WVector8_1(k)+WVector8_2(k));

        // 20

        WVector20(k) = getDivisionRecurrenceRelation(XMatrix.row(26),XMatrix.row(20),WVector20,k);


        UMatrix(20,k) = 0.5*WVector20(k);

        // 35

        WVector35_1(k) = getMultiplicationRecurrenceRelation(WVector4_6,XMatrix.row(25),k);
        WVector35_2(k) = getMultiplicationRecurrenceRelation(WVector4_4,XMatrix.row(24),k);
        WVector35_3(k) = getMultiplicationRecurrenceRelation(WVector35_2,XMatrix.row(12),k);

        UMatrix(35,k) = rotationalVelocity*(WVector35_1(k)-WVector35_3(k));

        // 9

        W9IntermediateVector(k) = getMultiplicationRecurrenceRelation(XMatrix.row(9),UMatrix.row(8),k);

        WVector9(k) = getDivisionRecurrenceRelation(W9IntermediateVector,XMatrix.row(8),WVector9,k);

        UMatrix(9,k) = 1.5*WVector9(k);

        // 21

        WVector21_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(4),UMatrix.row(4),k);
        WVector21_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(5),UMatrix.row(5),k);
        WVector21_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(6),UMatrix.row(6),k);


        UMatrix(21,k) = 2*(WVector21_1(k)+WVector21_2(k)+WVector21_3(k));

        // 26

        WVector26_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(1),UMatrix.row(4),k);
        WVector26_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(2),UMatrix.row(5),k);
        WVector26_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(3),UMatrix.row(6),k);

        UMatrix(26,k) = 2*XMatrix(21,k)+2*(WVector26_1(k)+WVector26_2(k)+WVector26_3(k));

        // 31

        UMatrix(31,k) = UMatrix(20,k);

        // 45

        WVector45_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(1),UMatrix.row(5),k);
        WVector45_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(2),UMatrix.row(4),k);
        WVector45_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(1),XMatrix.row(5),k);
        WVector45_4(k) = getMultiplicationRecurrenceRelation(XMatrix.row(2),XMatrix.row(4),k);
        WVector45_5(k) = getMultiplicationRecurrenceRelation(XMatrix.row(19),XMatrix.row(19),k);
        WVector45_6(k) = getMultiplicationRecurrenceRelation(XMatrix.row(19),(WVector45_1-WVector45_2),k);
        WVector45_7(k) = getMultiplicationRecurrenceRelation(UMatrix.row(19),(WVector45_3-WVector45_4),k);
        WVector45_8(k) = getDivisionRecurrenceRelation((WVector45_6-WVector45_7),WVector45_5,WVector45_8,k);


        UMatrix(45,k) = WVector45_8(k);

        // 25

        WVector25_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(8),UMatrix.row(26),k);
        WVector25_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(26),XMatrix.row(26),k);
        WVector25_3(k) = getDivisionRecurrenceRelation((2*WVector25_1-WVector25_2),XMatrix.row(9),WVector25_3,k);

        UMatrix(25,k) = 0.25*WVector25_3(k);

        // 30

        WVector30_1(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_1,9,k);
        WVector30_2(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_2,8,k);
        WVector30_3(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_3,7,k);
        WVector30_4(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_4,6,k);
        WVector30_5(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_5,5,k);
        WVector30_6(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_6,4,k);
        WVector30_7(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_7,3,k);
        WVector30_8(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_8,2,k);


        WVector30_9(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),(10*densityPolyCoefficients(10)*WVector30_1+9*densityPolyCoefficients(9)*WVector30_2+8*densityPolyCoefficients(8)*WVector30_3+
                                          7*densityPolyCoefficients(7)*WVector30_4+6*densityPolyCoefficients(6)*WVector30_5+5*densityPolyCoefficients(5)*WVector30_6+
                                          4*densityPolyCoefficients(4)*WVector30_7+3*densityPolyCoefficients(3)*WVector30_8+2*densityPolyCoefficients(2)*XMatrix.row(31)+densityPolyCoefficient_1),k);
                            // change to Eigen::VectorXd

        UMatrix(30,k) = WVector30_9(k);

        // 34
        // First it has to be determined which equation has to be used and which values have to be computed

        if (temperatureAltitudeRanges(0,0)<=XMatrix(31,0) && XMatrix(31,0) < temperatureAltitudeRanges(0,1)){

            UMatrix(34,k) = temperaturePolyCoefficients(0,1)*UMatrix(31,k);
        }

        else if (temperatureAltitudeRanges(1,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(1,1)){

            WVector34_2(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),(3*temperaturePolyCoefficients(1,2)*WVector30_8+2*temperaturePolyCoefficients(1,1)*XMatrix.row(31)+temperaturePolyCoefficient_1_2),k);
        // change to Eigen::VectorXd

            UMatrix(34,k) = WVector34_2(k);
        }

        else if (temperatureAltitudeRanges(2,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(2,1)){

            WVector34_3(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),(6*temperaturePolyCoefficients(2,5)*WVector30_5+5*temperaturePolyCoefficients(2,4)*WVector30_6+4*temperaturePolyCoefficients(2,3)*WVector30_7+
                                                                                  3*temperaturePolyCoefficients(2,2)*WVector30_8+2*temperaturePolyCoefficients(2,1)*XMatrix.row(31)+temperaturePolyCoefficient_1_3),k);
            // change to Eigen::VectorXd

            UMatrix(34,k) = WVector34_3(k);

            }

        else if (temperatureAltitudeRanges(3,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(3,1)){

            WVector34_4(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),(8*temperaturePolyCoefficients(3,7)*WVector30_3+7*temperaturePolyCoefficients(3,6)*WVector30_4+
                                                                     6*temperaturePolyCoefficients(3,5)*WVector30_5+5*temperaturePolyCoefficients(3,4)*WVector30_6+4*temperaturePolyCoefficients(3,3)*WVector30_7+
                                                                                  3*temperaturePolyCoefficients(3,2)*WVector30_8+2*temperaturePolyCoefficients(3,1)*XMatrix.row(31)+temperaturePolyCoefficient_1_4),k);
            // change to Eigen::VectorXd

            UMatrix(34,k) = WVector34_4(k);


        }
        else if (temperatureAltitudeRanges(4,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(4,1)){

            UMatrix(34,k) = 0;


    };




        // 36

        WVector36(k) = getDivisionRecurrenceRelation(UMatrix.row(21),XMatrix.row(36),WVector36,k);

        UMatrix(36,k) = 0.5*WVector36(k);

        // 46

        UMatrix(46,k) = UMatrix(45,k);

        // 47

        WVector47_1(k) =  getMultiplicationRecurrenceRelation(UMatrix.row(45),WVector4_6,k);
        WVector47_2(k) =  getMultiplicationRecurrenceRelation(UMatrix.row(12),WVector4_4,k);
        WVector47_3(k) =  getMultiplicationRecurrenceRelation(XMatrix.row(45),WVector47_2,k);

        UMatrix(47,k) = WVector47_1(k)-WVector47_3(k);

        // 24

        WVector24_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(6),UMatrix.row(20),k);
        WVector24_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(20),UMatrix.row(6),k);
        WVector24_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(25),XMatrix.row(6),k);
        WVector24_4(k) = getMultiplicationRecurrenceRelation(XMatrix.row(3),UMatrix.row(25),k);
        WVector24_5(k) = getDivisionRecurrenceRelation((WVector24_1+WVector24_2-WVector24_3-WVector24_4),WVector12_6,WVector24_5,k);
        WVector24_6(k) = getMultiplicationRecurrenceRelation(XMatrix.row(3),WVector24_4,k);
        WVector24_7(k) = getPowerRecurrenceRelation(XMatrix.row(20),WVector24_7,3,k);
        WVector24_8(k) = getDivisionRecurrenceRelation(WVector8_3,XMatrix.row(8),WVector24_8,k);
        WVector24_9(k) = getDivisionRecurrenceRelation(WVector24_6,WVector24_7,WVector24_9,k);
        WVector24_10(k) = getMultiplicationRecurrenceRelation((2*WVector24_9-2*WVector24_8),(WVector12_1-WVector12_2),k);

        WVector24_11(k) = getMultiplicationRecurrenceRelation(UMatrix.row(8),(WVector12_1-WVector12_2),k);
        WVector24_12(k) = getPowerRecurrenceRelation((onesVector-WVector12_4),WVector24_12,1.5,k);                                    // change to Eigen::VectorXd --- check
        WVector24_13(k) = getMultiplicationRecurrenceRelation(XMatrix.row(8),XMatrix.row(8),k);
        WVector24_14(k) = getMultiplicationRecurrenceRelation(XMatrix.row(8),WVector24_12,k);
        WVector24_15(k) = getMultiplicationRecurrenceRelation(WVector24_13,WVector12_5,k);
        WVector24_16(k) = getDivisionRecurrenceRelation(WVector24_10,WVector24_14,WVector24_16,k);
        WVector24_17(k) = getDivisionRecurrenceRelation(WVector24_11,WVector24_15,WVector24_17,k);

        UMatrix(24,k) = WVector24_5(k)-0.5*WVector24_16(k)-WVector24_17(k);


        // 28

        WVector28(k) = getMultiplicationRecurrenceRelation(UMatrix.row(30),XMatrix.row(28),k);

        UMatrix(28,k) = WVector28(k);

        // 33

        WVector33(k) = getDivisionRecurrenceRelation(UMatrix.row(34),XMatrix.row(33),WVector33,k);

        UMatrix(33,k) = WVector33(k)*(adiabeticIndex*specificGasConstant)/2;

        // 37

        WVector37_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(36),UMatrix.row(25),k);
        WVector37_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(25),UMatrix.row(36),k);
        WVector37_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(36),XMatrix.row(36),k);
        WVector37_4(k) = getDivisionRecurrenceRelation((WVector37_1-WVector37_2),WVector37_3,WVector37_4,k);

        UMatrix(37,k) = WVector37_4(k);

        // 41

        WVector40_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(36),UMatrix.row(35),k);
        WVector40_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(35),UMatrix.row(36),k);

        UMatrix(41,k) = WVector41_1(k)+WVector41_2(k);

        // 48

        WVector48_1(k) =  getMultiplicationRecurrenceRelation(UMatrix.row(46),WVector4_6,k);
        WVector48_2(k) =  getMultiplicationRecurrenceRelation(XMatrix.row(46),WVector47_2,k);

        UMatrix(48,k) = WVector48_1(k)-WVector48_2(k);

        // 13

        WVector13_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(24),UMatrix.row(48),k);
        WVector13_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(48),UMatrix.row(24),k);
        WVector13_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(48),XMatrix.row(48),k);
        WVector13_4(k) = getMultiplicationRecurrenceRelation(XMatrix.row(24),XMatrix.row(24),k);
        WVector13_5(k) = getDivisionRecurrenceRelation((WVector13_1-WVector13_2),(WVector13_3+WVector13_4),WVector13_5,k);


        UMatrix(13,k) = WVector13_5(k);

        // 38

        WVector38_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(37),XMatrix.row(37),k);
        WVector38_2(k) = getPowerRecurrenceRelation((onesVector-WVector38_1),WVector38_2,0.5,k);                         // change to Eigen::VectorXd --- check
        WVector38_3(k) = getDivisionRecurrenceRelation(UMatrix.row(37),WVector38_2,WVector38_3,k);


        UMatrix(38,k) = WVector38_3(k);

        // 40

        WVector40_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(24),UMatrix.row(47),k);
        WVector40_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(47),UMatrix.row(24),k);
        WVector40_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(47),XMatrix.row(47),k);
        WVector40_4(k) = getDivisionRecurrenceRelation((WVector40_1-WVector40_2),(WVector40_3+WVector13_4),WVector40_4,k);



        UMatrix(40,k) = WVector40_4(k);

        // 42

        WVector42_1(k) = getCosineRecurrenceRelation(XMatrix.row(38),WVector42_3,k);
        WVector42_2(k) = getCosineRecurrenceRelation(XMatrix.row(40),WVector42_4,k);
        WVector42_3(k) = getSineRecurrenceRelation(XMatrix.row(38),WVector42_1,k);
        WVector42_4(k) = getSineRecurrenceRelation(XMatrix.row(40),WVector42_2,k);
        WVector42_5(k) = getMultiplicationRecurrenceRelation(WVector42_2,UMatrix.row(40),k);
        WVector42_6(k) = getMultiplicationRecurrenceRelation(WVector42_3,UMatrix.row(38),k);
        WVector42_7(k) = getMultiplicationRecurrenceRelation(WVector42_1,WVector42_5,k);
        WVector42_8(k) = getMultiplicationRecurrenceRelation(WVector42_6,WVector42_4,k);

        UMatrix(42,k) = WVector42_7(k) - WVector42_8(k);

        // 43

        WVector43_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(41),UMatrix.row(42),k);
        WVector43_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(42),UMatrix.row(41),k);

        UMatrix(43,k) = WVector43_1(k) + WVector43_2(k);

        // 15

        WVector15_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(35),UMatrix.row(35),k);
        WVector15_2(k) = getDivisionRecurrenceRelation((2*WVector15_1+UMatrix.row(21)-2*UMatrix.row(43)),XMatrix.row(15),WVector15_2,k);

        UMatrix(15,k) = 0.5*WVector15_2(k);

        // 23

        WVector23_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(15),UMatrix.row(25),k);
        WVector23_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(25),UMatrix.row(15),k);
        WVector23_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(15),XMatrix.row(15),k);
        WVector23_4(k) = getDivisionRecurrenceRelation((WVector23_1-WVector23_2),WVector23_3,WVector23_4,k);


        UMatrix(23,k) = WVector23_4(k);

        // 32

        WVector32_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(33),UMatrix.row(15),k);
        WVector32_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(15),UMatrix.row(33),k);
        WVector32_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(33),XMatrix.row(33),k);
        WVector32_4(k) = getDivisionRecurrenceRelation((WVector32_1-WVector32_2),WVector32_3,WVector32_4,k);

        UMatrix(32,k) = WVector32_4(k);

        // 14

        WVector14_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(23),XMatrix.row(23),k);
        WVector14_2(k) = getPowerRecurrenceRelation((onesVector-WVector14_1),WVector14_2,0.5,k);                         // change to Eigen::VectorXd --- check
        WVector14_3(k) = getDivisionRecurrenceRelation(UMatrix.row(23),WVector14_2,WVector14_3,k);


        UMatrix(14,k) = WVector14_3(k);

        // 29

        for (int i=0; i < 5+1; i++){

            if (dragCoefficientMachRanges(i,0) <= XMatrix(32,0) && XMatrix(32,0) < dragCoefficientMachRanges(i,1)){

               sectionCD = i;
            }

        };

        UMatrix(29,k) = dragCoefficientPolyCoefficients(sectionCD,1)*UMatrix(32,k);


        // 16

        WVector16(k) = getMultiplicationRecurrenceRelation(WVector4_7,UMatrix.row(14),k);

        UMatrix(16,k) = -WVector16(k);

        // 27

        WVector27_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(29),UMatrix.row(28),k);
        WVector27_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(28),UMatrix.row(29),k);
        WVector27_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(28),XMatrix.row(29),k);
        WVector27_4(k) = getMultiplicationRecurrenceRelation(XMatrix.row(15),(WVector27_1+WVector27_2),k);
        WVector27_5(k) = getMultiplicationRecurrenceRelation(WVector27_3,UMatrix.row(15),k);
        WVector27_6(k) = getMultiplicationRecurrenceRelation(XMatrix.row(15),(WVector27_4+WVector27_5),k);

        UMatrix(27,k) = 0.5*referenceArea*WVector27_6(k);



        /// Compute all auxiliary Equation coefficients from the Derivative coefficients ///

        for (int i = 1; i< initialEquationsVector_.size(); i++){

            XMatrix(i,k+1) = UMatrix(i,k)/(k+1);
        }


//        std::cout<<"Number of runs "<<count<<std::endl;

        count++;

//        std::cout<<"So far so good 5"<<std::endl;

//        std::cout<<"Length of W = "<<WVector4_1.size()<<std::endl;
//        std::cout<<"Length of U = "<<UMatrix.row(0).size()<<std::endl;

} // end of for loop with k till K



//    std::cout<<"So far so good 6"<<std::endl;

    /// Set return matrix ///

    Eigen::MatrixXd stateTaylorCoefficients = Eigen::MatrixXd::Zero(8,maxOrder+1);

    stateTaylorCoefficients.row(1) = XMatrix.row(1);
    stateTaylorCoefficients.row(2) = XMatrix.row(2);
    stateTaylorCoefficients.row(3) = XMatrix.row(3);
    stateTaylorCoefficients.row(4) = XMatrix.row(4);
    stateTaylorCoefficients.row(5) = XMatrix.row(5);
    stateTaylorCoefficients.row(6) = XMatrix.row(6);
    stateTaylorCoefficients.row(7) = XMatrix.row(7);


    /// Debug ///
/*
//    std::cout<<"So far so good 7"<<std::endl;

    std::cout<<"stateTaylorCoefficients = "<<stateTaylorCoefficients<<std::endl;

//    Eigen::VectorXd debugVecStateTaylorCoefficients = stateTaylorCoefficients.row(1);

    std::cout<<"Length of stateTaylorCoefficients.row(1) = "<<stateTaylorCoefficients.row(1).size()<<std::endl;
    std::cout<<"Length of XMatrix.row(1) = "<<XMatrix.row(1).size()<<std::endl;


    // So it works until it has to be returned... So... delete everything? :S
//*/
    return stateTaylorCoefficients;

} // end of function




