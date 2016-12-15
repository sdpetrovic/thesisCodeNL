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
 *      160518    S.D. Petrovic     Fixed the mistake I made with the transformation matrix T_IB which resulted in a mistake in u6
 *      160520    S.D. Petrovic     Fixed the mistake where I had written W40_1 and W40_2 instead of W41_1 and W41_2
 *      160526    S.D. Petrovic     Added W4,0 to be able to properly evaluate the recurrence relation of W4,2
 *      160602    S.D. Petrovic     Added the thrust recurrence relations
 *      160623    S.D. Petrovic     Corrected U15 and added corresponding auxiliary functions and recurrence relations
 *      160628    S.D. Petrovic     Rewrote all transformation angle equations!!
 *      160712    S.D. Petrovic     Tested new way of writing recurrence relations using Bergsma's thesis
 *
 *    References
 *
 *    Notes
 *
 */



#include "allRecurrenceRelationsSpherical.h"
#include "basicRecurrenceRelations.h"



Eigen::MatrixXd getTaylorCoefficients(const double adiabeticIndex_, const double specificGasConstant_, const double standardGravitationalParameter_, const double rotationalVelocity_, const double primeMeridianAngle_,
                      const double inertialFrameTime_, const double bodyReferenceRadius_, const Eigen::MatrixXd temperaturePolyCoefficients_, const Eigen::MatrixXd temperatureAltitudeRanges_,
                      const Eigen::VectorXd densityPolyCoefficients_, const double Thrust_, const double specificImpulse_,
                      const double referenceArea_, const Eigen::MatrixXd dragCoefficientPolyCoefficients_, const Eigen::MatrixXd dragCoefficientMachRanges_,
        const Eigen::VectorXd thrustAccelerationsBframe_,
        const Eigen::MatrixXd thrustAzimuthMatrix_,
        const Eigen::MatrixXd thrustElevationMatrix_,
        const Eigen::VectorXd initialEquationsVector_,
        const Eigen::VectorXd initialDerivativesVector_,
        const Eigen::MatrixXd initialFunctionsMatrix_,
        const double currentTime_,
        const int maxOrder_){


    // Create new variables to use in this function and be deleted afterwards... I hope... Update: Nope... Update2: Yep :)
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
    const Eigen::MatrixXd thrustAzimuthMatrix = thrustAzimuthMatrix_;
    const Eigen::MatrixXd thrustElevationMatrix = thrustElevationMatrix_;
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

    for (int i = 1; i < initialEquationsVector.size(); i++){                      // Fill the matrices. Please note that the number allocation and position allocation in the matrix/vector should always be the same.
                                                                                  // This means that x1 will be in position (1) of the vector (which is the second vector entry!)
      XMatrix(i,0) = initialEquationsVector(i);
      XMatrix(i,1) = initialDerivativesVector(i);
      UMatrix(i,0) = initialDerivativesVector(i);




    };


    // Create the auxiliary functions vectors

    Eigen::VectorXd WVector4_0 = Eigen::VectorXd::Zero(maxOrder);     // W4,0   // Added because of the mistake found in the recurrence relation of W4,2
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
    Eigen::VectorXd WVector4_25 = Eigen::VectorXd::Zero(maxOrder);     // W4,25
    Eigen::VectorXd WVector4_26 = Eigen::VectorXd::Zero(maxOrder);     // W4,26
    Eigen::VectorXd WVector4_27 = Eigen::VectorXd::Zero(maxOrder);     // W4,27
    Eigen::VectorXd WVector4_28 = Eigen::VectorXd::Zero(maxOrder);     // W4,28
    Eigen::VectorXd WVector4_29 = Eigen::VectorXd::Zero(maxOrder);     // W4,29
    Eigen::VectorXd WVector4_30 = Eigen::VectorXd::Zero(maxOrder);     // W4,30

    Eigen::VectorXd WVector4_31 = Eigen::VectorXd::Zero(maxOrder);     // W4,31
    Eigen::VectorXd WVector4_32 = Eigen::VectorXd::Zero(maxOrder);     // W4,32
    Eigen::VectorXd WVector4_33 = Eigen::VectorXd::Zero(maxOrder);     // W4,33
    Eigen::VectorXd WVector4_34 = Eigen::VectorXd::Zero(maxOrder);     // W4,34
    Eigen::VectorXd WVector4_35 = Eigen::VectorXd::Zero(maxOrder);     // W4,35
    Eigen::VectorXd WVector4_36 = Eigen::VectorXd::Zero(maxOrder);     // W4,36
    Eigen::VectorXd WVector4_37 = Eigen::VectorXd::Zero(maxOrder);     // W4,37
    Eigen::VectorXd WVector4_38 = Eigen::VectorXd::Zero(maxOrder);     // W4,38


    Eigen::VectorXd WVector5_1 = Eigen::VectorXd::Zero(maxOrder);     // W5,1
    Eigen::VectorXd WVector5_2 = Eigen::VectorXd::Zero(maxOrder);     // W5,2
    Eigen::VectorXd WVector5_3 = Eigen::VectorXd::Zero(maxOrder);     // W5,3
    Eigen::VectorXd WVector5_4 = Eigen::VectorXd::Zero(maxOrder);     // W5,4
    Eigen::VectorXd WVector5_5 = Eigen::VectorXd::Zero(maxOrder);     // W5,5
    Eigen::VectorXd WVector5_6 = Eigen::VectorXd::Zero(maxOrder);     // W5,6
    Eigen::VectorXd WVector5_7 = Eigen::VectorXd::Zero(maxOrder);     // W5,7
    Eigen::VectorXd WVector5_8 = Eigen::VectorXd::Zero(maxOrder);     // W5,8
    Eigen::VectorXd WVector5_9 = Eigen::VectorXd::Zero(maxOrder);     // W5,9
    Eigen::VectorXd WVector5_10 = Eigen::VectorXd::Zero(maxOrder);     // W5,10


    Eigen::VectorXd WVector6_0 = Eigen::VectorXd::Zero(maxOrder);     // W6,0  // Added because of the mistake found in the complete transformation matrix
    Eigen::VectorXd WVector6_1 = Eigen::VectorXd::Zero(maxOrder);     // W6,1
    Eigen::VectorXd WVector6_2 = Eigen::VectorXd::Zero(maxOrder);     // W6,2
    Eigen::VectorXd WVector6_3 = Eigen::VectorXd::Zero(maxOrder);     // W6,3
    Eigen::VectorXd WVector6_4 = Eigen::VectorXd::Zero(maxOrder);     // W6,4
    Eigen::VectorXd WVector6_5 = Eigen::VectorXd::Zero(maxOrder);     // W6,5
    Eigen::VectorXd WVector6_6 = Eigen::VectorXd::Zero(maxOrder);     // W6,6
    Eigen::VectorXd WVector6_7 = Eigen::VectorXd::Zero(maxOrder);     // W6,7
    Eigen::VectorXd WVector6_8 = Eigen::VectorXd::Zero(maxOrder);     // W6,8


    Eigen::VectorXd WVector8_1 = Eigen::VectorXd::Zero(maxOrder);     // W8,1
    Eigen::VectorXd WVector8_2 = Eigen::VectorXd::Zero(maxOrder);     // W8,2
    Eigen::VectorXd WVector8_3 = Eigen::VectorXd::Zero(maxOrder);     // W8,3


    Eigen::VectorXd WVector9 = Eigen::VectorXd::Zero(maxOrder);     // W9


    Eigen::VectorXd WVector11_0 = Eigen::VectorXd::Zero(maxOrder);     // W11,0
    Eigen::VectorXd WVector11_1 = Eigen::VectorXd::Zero(maxOrder);     // W11,1
    Eigen::VectorXd WVector11_2 = Eigen::VectorXd::Zero(maxOrder);     // W11,2
    Eigen::VectorXd WVector11_3 = Eigen::VectorXd::Zero(maxOrder);     // W11,3

    Eigen::VectorXd WVector12_1 = Eigen::VectorXd::Zero(maxOrder);     // W12,1
//    Eigen::VectorXd WVector12_2 = Eigen::VectorXd::Zero(maxOrder);     // W12,2

    Eigen::VectorXd WVector13_0 = Eigen::VectorXd::Zero(maxOrder);     // W13,0
    Eigen::VectorXd WVector13_1 = Eigen::VectorXd::Zero(maxOrder);     // W13,1
    Eigen::VectorXd WVector13_2 = Eigen::VectorXd::Zero(maxOrder);     // W13,2
    Eigen::VectorXd WVector13_3 = Eigen::VectorXd::Zero(maxOrder);     // W13,3
    Eigen::VectorXd WVector13_4 = Eigen::VectorXd::Zero(maxOrder);     // W13,4
    Eigen::VectorXd WVector13_5 = Eigen::VectorXd::Zero(maxOrder);     // W13,5
    Eigen::VectorXd WVector13_6 = Eigen::VectorXd::Zero(maxOrder);     // W13,6
    Eigen::VectorXd WVector13_7 = Eigen::VectorXd::Zero(maxOrder);     // W13,7
    Eigen::VectorXd WVector13_8 = Eigen::VectorXd::Zero(maxOrder);     // W13,8
    Eigen::VectorXd WVector13_9 = Eigen::VectorXd::Zero(maxOrder);     // W13,9
//    Eigen::VectorXd WVector13_10 = Eigen::VectorXd::Zero(maxOrder);     // W13,10
//    Eigen::VectorXd WVector13_11 = Eigen::VectorXd::Zero(maxOrder);     // W13,11
//    Eigen::VectorXd WVector13_12 = Eigen::VectorXd::Zero(maxOrder);     // W13,12

    Eigen::VectorXd WVector14_0 = Eigen::VectorXd::Zero(maxOrder);     // W14,0
    Eigen::VectorXd WVector14_1 = Eigen::VectorXd::Zero(maxOrder);     // W14,1
    Eigen::VectorXd WVector14_2 = Eigen::VectorXd::Zero(maxOrder);     // W14,2
    Eigen::VectorXd WVector14_3 = Eigen::VectorXd::Zero(maxOrder);     // W14,3
    Eigen::VectorXd WVector14_4 = Eigen::VectorXd::Zero(maxOrder);     // W14,4
    Eigen::VectorXd WVector14_5 = Eigen::VectorXd::Zero(maxOrder);     // W14,5
    Eigen::VectorXd WVector14_6 = Eigen::VectorXd::Zero(maxOrder);     // W14,6
    Eigen::VectorXd WVector14_7 = Eigen::VectorXd::Zero(maxOrder);     // W14,7
//    Eigen::VectorXd WVector14_8 = Eigen::VectorXd::Zero(maxOrder);     // W14,8
//    Eigen::VectorXd WVector14_9 = Eigen::VectorXd::Zero(maxOrder);     // W14,9


    Eigen::VectorXd WVector15_0 = Eigen::VectorXd::Zero(maxOrder);     // W15,0
    Eigen::VectorXd WVector15_1 = Eigen::VectorXd::Zero(maxOrder);     // W15,1
    Eigen::VectorXd WVector15_2 = Eigen::VectorXd::Zero(maxOrder);     // W15,2
    Eigen::VectorXd WVector15_3 = Eigen::VectorXd::Zero(maxOrder);     // W15,3
    Eigen::VectorXd WVector15_4 = Eigen::VectorXd::Zero(maxOrder);     // W15,4
    Eigen::VectorXd WVector15_5 = Eigen::VectorXd::Zero(maxOrder);     // W15,5
    Eigen::VectorXd WVector15_6 = Eigen::VectorXd::Zero(maxOrder);     // W15,6

    Eigen::VectorXd WVector16_1 = Eigen::VectorXd::Zero(maxOrder);     // W16,1


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
    Eigen::VectorXd WVector30_19 = Eigen::VectorXd::Zero(maxOrder);     // W30,19
    Eigen::VectorXd WVector30_9 = Eigen::VectorXd::Zero(maxOrder);     // W30,9


    Eigen::VectorXd WVector32_1 = Eigen::VectorXd::Zero(maxOrder);     // W32,1
    Eigen::VectorXd WVector32_2 = Eigen::VectorXd::Zero(maxOrder);     // W32,2
    Eigen::VectorXd WVector32_3 = Eigen::VectorXd::Zero(maxOrder);     // W32,3
    Eigen::VectorXd WVector32_4 = Eigen::VectorXd::Zero(maxOrder);     // W32,4


    Eigen::VectorXd WVector33 = Eigen::VectorXd::Zero(maxOrder);     // W33


    Eigen::VectorXd WVector34_12 = Eigen::VectorXd::Zero(maxOrder);     // W34,12
    Eigen::VectorXd WVector34_2 = Eigen::VectorXd::Zero(maxOrder);     // W34,2
    Eigen::VectorXd WVector34_13 = Eigen::VectorXd::Zero(maxOrder);     // W34,13
    Eigen::VectorXd WVector34_3 = Eigen::VectorXd::Zero(maxOrder);     // W34,3
    Eigen::VectorXd WVector34_14 = Eigen::VectorXd::Zero(maxOrder);     // W34,14
    Eigen::VectorXd WVector34_4 = Eigen::VectorXd::Zero(maxOrder);     // W34,4




    // And fill them

    WVector4_0(0) = initialFunctionsMatrix(4,0);     // W4,0   // Added because of the mistake found in the recurrence relation of W4,2
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
    WVector4_25(0) = initialFunctionsMatrix(4,25);     // W4,25
    WVector4_26(0) = initialFunctionsMatrix(4,26);     // W4,26
    WVector4_27(0) = initialFunctionsMatrix(4,27);     // W4,27
    WVector4_28(0) = initialFunctionsMatrix(4,28);     // W4,28
    WVector4_29(0) = initialFunctionsMatrix(4,29);     // W4,29
    WVector4_30(0) = initialFunctionsMatrix(4,30);     // W4,30

    WVector4_31(0) = initialFunctionsMatrix(4,31);     // W4,31
    WVector4_32(0) = initialFunctionsMatrix(4,32);     // W4,32
    WVector4_33(0) = initialFunctionsMatrix(4,33);     // W4,33
    WVector4_34(0) = initialFunctionsMatrix(4,34);     // W4,34
    WVector4_35(0) = initialFunctionsMatrix(4,35);     // W4,35
    WVector4_36(0) = initialFunctionsMatrix(4,36);     // W4,36
    WVector4_37(0) = initialFunctionsMatrix(4,37);     // W4,37
    WVector4_38(0) = initialFunctionsMatrix(4,38);     // W4,38


    WVector5_1(0) = initialFunctionsMatrix(5,1);     // W5,1
    WVector5_2(0) = initialFunctionsMatrix(5,2);     // W5,2
    WVector5_3(0) = initialFunctionsMatrix(5,3);     // W5,3
    WVector5_4(0) = initialFunctionsMatrix(5,4);     // W5,4
    WVector5_5(0) = initialFunctionsMatrix(5,5);     // W5,5
    WVector5_6(0) = initialFunctionsMatrix(5,6);     // W5,6
    WVector5_7(0) = initialFunctionsMatrix(5,7);     // W5,7
    WVector5_8(0) = initialFunctionsMatrix(5,8);     // W5,8
    WVector5_9(0) = initialFunctionsMatrix(5,9);     // W5,9
    WVector5_10(0) = initialFunctionsMatrix(5,10);     // W5,10


    WVector6_0(0) = initialFunctionsMatrix(6,0);     // W6,0  // Added because of the mistake found in the complete transformation matrix
    WVector6_1(0) = initialFunctionsMatrix(6,1);     // W6,1
    WVector6_2(0) = initialFunctionsMatrix(6,2);     // W6,2
    WVector6_3(0) = initialFunctionsMatrix(6,3);     // W6,3
    WVector6_4(0) = initialFunctionsMatrix(6,4);     // W6,4
    WVector6_5(0) = initialFunctionsMatrix(6,5);     // W6,5
    WVector6_6(0) = initialFunctionsMatrix(6,6);     // W6,6
    WVector6_7(0) = initialFunctionsMatrix(6,7);     // W6,7
    WVector6_8(0) = initialFunctionsMatrix(6,8);     // W6,8


    WVector8_1(0) = initialFunctionsMatrix(8,1);     // W8,1
    WVector8_2(0) = initialFunctionsMatrix(8,2);     // W8,2
    WVector8_3(0) = initialFunctionsMatrix(8,3);     // W8,3


    WVector9(0) = initialFunctionsMatrix(9,1);     // W9


    WVector11_0(0) = initialFunctionsMatrix(11,0);     // W11,0
    WVector11_1(0) = initialFunctionsMatrix(11,1);     // W11,1
    WVector11_2(0) = initialFunctionsMatrix(11,2);     // W11,2
    WVector11_3(0) = initialFunctionsMatrix(11,3);     // W11,3

    WVector12_1(0) = initialFunctionsMatrix(12,1);     // W12,1
//    WVector12_2(0) = initialFunctionsMatrix(12,2);     // W12,2

    WVector13_0(0) = initialFunctionsMatrix(13,0);     // W13,0
    WVector13_1(0) = initialFunctionsMatrix(13,1);     // W13,1
    WVector13_2(0) = initialFunctionsMatrix(13,2);     // W13,2
    WVector13_3(0) = initialFunctionsMatrix(13,3);     // W13,3
    WVector13_4(0) = initialFunctionsMatrix(13,4);     // W13,4
    WVector13_5(0) = initialFunctionsMatrix(13,5);     // W13,5
    WVector13_6(0) = initialFunctionsMatrix(13,6);     // W13,6
    WVector13_7(0) = initialFunctionsMatrix(13,7);     // W13,7
    WVector13_8(0) = initialFunctionsMatrix(13,8);     // W13,8
    WVector13_9(0) = initialFunctionsMatrix(13,9);     // W13,9
//    WVector13_10(0) = initialFunctionsMatrix(13,10);     // W13,10
//    WVector13_11(0) = initialFunctionsMatrix(13,11);     // W13,11
//    WVector13_12(0) = initialFunctionsMatrix(13,12);     // W13,12

    WVector14_0(0) = initialFunctionsMatrix(14,0);     // W14,1
    WVector14_1(0) = initialFunctionsMatrix(14,1);     // W14,1
    WVector14_2(0) = initialFunctionsMatrix(14,2);     // W14,2
    WVector14_3(0) = initialFunctionsMatrix(14,3);     // W14,3
    WVector14_4(0) = initialFunctionsMatrix(14,4);     // W14,4
    WVector14_5(0) = initialFunctionsMatrix(14,5);     // W14,5
    WVector14_6(0) = initialFunctionsMatrix(14,6);     // W14,6
    WVector14_7(0) = initialFunctionsMatrix(14,7);     // W14,7
//    WVector14_8(0) = initialFunctionsMatrix(14,8);     // W14,8
//    WVector14_9(0) = initialFunctionsMatrix(14,9);     // W14,9

    WVector15_0(0) = initialFunctionsMatrix(15,0);     // W15,0
    WVector15_1(0) = initialFunctionsMatrix(15,1);     // W15,1
    WVector15_2(0) = initialFunctionsMatrix(15,2);     // W15,2
    WVector15_3(0) = initialFunctionsMatrix(15,3);     // W15,3
    WVector15_4(0) = initialFunctionsMatrix(15,4);     // W15,4
    WVector15_5(0) = initialFunctionsMatrix(15,5);     // W15,5
    WVector15_6(0) = initialFunctionsMatrix(15,6);     // W15,6

    WVector16_1(0) = initialFunctionsMatrix(16,1);     // W16,1


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
    WVector30_19(0) = initialFunctionsMatrix(30,19);     // W30,19
    WVector30_9(0) = initialFunctionsMatrix(30,9);     // W30,9


    WVector32_1(0) = initialFunctionsMatrix(32,1);     // W32,1
    WVector32_2(0) = initialFunctionsMatrix(32,2);     // W32,2
    WVector32_3(0) = initialFunctionsMatrix(32,3);     // W32,3
    WVector32_4(0) = initialFunctionsMatrix(32,4);     // W32,4


    WVector33(0) = initialFunctionsMatrix(33,1);     // W33


    WVector34_12(0) = initialFunctionsMatrix(34,12);     // W34,12
    WVector34_2(0) = initialFunctionsMatrix(34,2);     // W34,2
    WVector34_13(0) = initialFunctionsMatrix(34,13);     // W34,13
    WVector34_3(0) = initialFunctionsMatrix(34,3);     // W34,3
    WVector34_14(0) = initialFunctionsMatrix(34,14);     // W34,14
    WVector34_4(0) = initialFunctionsMatrix(34,4);     // W34,4

    /// Debug ///
    Eigen::VectorXd WVectorTest = Eigen::VectorXd::Zero(maxOrder);     // Wtest
    WVectorTest(0) = XMatrix(15,0)*WVector4_38(0)/XMatrix(16,0);
    /// Debug ///



    ////// Declaring other Vectors and variables ///////


    Eigen::VectorXd W9IntermediateVector = Eigen::VectorXd::Zero(maxOrder);         // Intermediate vector created to be able to use the basic recurrence relations (should have done this from the start...)
    W9IntermediateVector(0) = XMatrix(9,0)*UMatrix(8,0);            // And setting the first value

//    Eigen::VectorXd W4_25IntermediateVector = Eigen::VectorXd::Zero(maxOrder);      // Intermediate vector created for the recurrence relation of thrust
//    W4_25IntermediateVector(0) = 1/XMatrix(7,0);

//    std::cout<<"W4_25IntermediateVector = "<<W4_25IntermediateVector<<std::endl;


    Eigen::VectorXd onesVector = Eigen::VectorXd::Zero(maxOrder);       // Vector used to represent the int 1 in cases where the entire vector has to be substracted from 1(for instance)

    Eigen::VectorXd densityPolyCoefficient_1 = Eigen::VectorXd::Zero(maxOrder);             // Vectors used to represent the polynomials when required for addition of vectors
    Eigen::VectorXd temperaturePolyCoefficient_1_1 = Eigen::VectorXd::Zero(maxOrder);
    Eigen::VectorXd temperaturePolyCoefficient_2_1 = Eigen::VectorXd::Zero(maxOrder);
    Eigen::VectorXd temperaturePolyCoefficient_3_1 = Eigen::VectorXd::Zero(maxOrder);
    Eigen::VectorXd ThrustVector = Eigen::VectorXd::Zero(maxOrder);
    Eigen::VectorXd thrustAzimuthVector = Eigen::VectorXd::Zero(maxOrder);
    Eigen::VectorXd thrustElevationVector = Eigen::VectorXd::Zero(maxOrder);

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

//            onesVector(i) = 1;
            densityPolyCoefficient_1(i) = densityPolyCoefficients(1);
            temperaturePolyCoefficient_1_1(i) = temperaturePolyCoefficients(1,1);
            temperaturePolyCoefficient_2_1(i) = temperaturePolyCoefficients(2,1);
            temperaturePolyCoefficient_3_1(i) = temperaturePolyCoefficients(3,1);
            ThrustVector(i) = Thrust_;
            if (i==0){
            thrustAzimuthVector(i) = thrustAzimuthMatrix(0,2);
            thrustElevationVector(i) = thrustElevationMatrix(0,2);

//            onesVector(i) = 1;
};

        }


//    std::cout<<"So far so good 3"<<std::endl;

//        // 1
//        UMatrix(1,k) = XMatrix(4,k);

//        // 2
//        UMatrix(2,k) = XMatrix(5,k);

//        // 3
//        UMatrix(3,k) = XMatrix(6,k);

        // 4

//        WVector4_0(k) = getDivisionRecurrenceRelation(XMatrix.row(27),XMatrix.row(7),WVector4_0,k); // Added because of the mistake found in the recurrence relation of W4,2
//        WVector4_1(k) = getDivisionRecurrenceRelation(XMatrix.row(1),XMatrix.row(9),WVector4_1,k); //  XMatrix(1,k)/XMatrix(9,k);
//        WVector4_2(k) = thrustAccelerationsBframe(0)-WVector4_0(k);                                     /// This is wrong! For the division only the previous division should be taken. So an extra WVector should be added!!  Update: done

        // Thrust additions
        WVector4_37(k) = getPowerRecurrenceRelation(XMatrix.row(7),WVector4_37,-1.0,k);
//        WVector4_25(k) = getDivisionRecurrenceRelation(ThrustVector,XMatrix.row(7),WVector4_25,k);
        WVector4_25(k) = Thrust_*WVector4_37(k);
        WVector4_26(k) = getCosineRecurrenceRelation(thrustAzimuthVector,WVector4_28,k);
        WVector4_27(k) = getCosineRecurrenceRelation(thrustElevationVector,WVector4_29,k);
        WVector4_28(k) = getSineRecurrenceRelation(thrustAzimuthVector,WVector4_26,k);
        WVector4_29(k) = getSineRecurrenceRelation(thrustElevationVector,WVector4_27,k);
        WVector4_30(k) = getMultiplicationRecurrenceRelation(WVector4_26,WVector4_27,k);
        WVector4_31(k) = getMultiplicationRecurrenceRelation(WVector4_28,WVector4_27,k);
        WVector4_32(k) = getMultiplicationRecurrenceRelation(WVector4_25,WVector4_30,k);

////        WVector4_2(k) = WVector4_32(k)-WVector4_0(k);
////        WVector4_3(k) = getCosineRecurrenceRelation((XMatrix.row(10)+XMatrix.row(11)),WVector4_8,k);
//        WVector4_3(k) = getCosineRecurrenceRelation((XMatrix.row(49)),WVector4_8,k);
        WVector4_4(k) = getSineRecurrenceRelation(XMatrix.row(12),WVector4_6,k);
        WVector4_5(k) = getCosineRecurrenceRelation(XMatrix.row(13),WVector4_9,k);
        WVector4_6(k) = getCosineRecurrenceRelation(XMatrix.row(12),WVector4_4,k);
////        WVector4_7(k) = getSineRecurrenceRelation(XMatrix.row(14),XMatrix.row(16),k);
        WVector4_7(k) = getSineRecurrenceRelation(XMatrix.row(14),WVector4_38,k);
        WVector4_38(k) = getCosineRecurrenceRelation(XMatrix.row(14),WVector4_7,k);
////        WVector4_8(k) = getSineRecurrenceRelation((XMatrix.row(10)+XMatrix.row(11)),WVector4_3,k);
        WVector4_8(k) = getSineRecurrenceRelation((XMatrix.row(49)),WVector4_3,k);
        WVector4_9(k) = getSineRecurrenceRelation(XMatrix.row(13),WVector4_5,k);
//        WVector4_10(k) = getMultiplicationRecurrenceRelation(WVector4_4,WVector4_5,k);
//        WVector4_11(k) = getMultiplicationRecurrenceRelation(WVector4_6,WVector4_7,k);
////        WVector4_12(k) = getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(16),k);
//        WVector4_12(k) = getMultiplicationRecurrenceRelation(WVector4_9,WVector4_38,k);
//        WVector4_13(k) = getMultiplicationRecurrenceRelation(WVector4_4,WVector4_9,k);
//        WVector4_14(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_5,k);
////        WVector4_15(k) = getMultiplicationRecurrenceRelation(WVector4_6,XMatrix.row(16),k);
//        WVector4_15(k) = getMultiplicationRecurrenceRelation(WVector4_6,WVector4_38,k);
//        WVector4_16(k) = getMultiplicationRecurrenceRelation(WVector4_9,WVector4_7,k);
////        WVector4_17(k) = getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(16),k);
//        WVector4_17(k) = getMultiplicationRecurrenceRelation(WVector4_10,WVector4_38,k);
//        WVector4_18(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_12,k);
//        WVector4_19(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_13,k);
//        WVector4_20(k) = getMultiplicationRecurrenceRelation(WVector4_10,WVector4_7,k);
//        WVector4_21(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_16,k);
//        WVector4_22(k) = getMultiplicationRecurrenceRelation(WVector4_3,(WVector4_11-WVector4_17),k);
//        WVector4_23(k) = getMultiplicationRecurrenceRelation(WVector4_3,(-WVector4_20-WVector4_15),k);
//        WVector4_24(k) = getMultiplicationRecurrenceRelation(WVector4_2,(WVector4_22-WVector4_18),k);

//        WVector4_33(k) = getMultiplicationRecurrenceRelation(WVector4_25,WVector4_31,k);
//        WVector4_34(k) = getMultiplicationRecurrenceRelation(WVector4_25,WVector4_29,k);
//        WVector4_35(k) = getMultiplicationRecurrenceRelation(WVector4_33,(WVector4_19-WVector4_14),k);
//        WVector4_36(k) = getMultiplicationRecurrenceRelation(WVector4_34,(WVector4_23-WVector4_21),k);

//        UMatrix(4,k) = -standardGravitationalParameter*WVector4_1(k)+WVector4_24(k)+WVector4_35(k)-WVector4_36(k);



//        // 5

////        WVector5_1(k) = getDivisionRecurrenceRelation(XMatrix.row(2),XMatrix.row(9),WVector5_1,k);
//        WVector5_2(k) = getMultiplicationRecurrenceRelation(WVector4_8,(WVector4_11-WVector4_17),k);
//        WVector5_3(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_12,k);
//        WVector5_4(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_13,k);
//        WVector5_5(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_5,k);
//        WVector5_6(k) = getMultiplicationRecurrenceRelation(WVector4_8,(-WVector4_20-WVector4_11),k);
//        WVector5_7(k) = getMultiplicationRecurrenceRelation(WVector4_3,WVector4_16,k);
//        WVector5_8(k) = getMultiplicationRecurrenceRelation(WVector4_2,(WVector5_2+WVector5_3),k);
//        WVector5_9(k) = getMultiplicationRecurrenceRelation(WVector4_33,(WVector5_4+WVector5_5),k);
//        WVector5_10(k) = getMultiplicationRecurrenceRelation(WVector4_34,(WVector5_6+WVector5_7),k);


////        UMatrix(5,k) = -standardGravitationalParameter*WVector5_1(k)+WVector5_8(k)+WVector5_9(k)-WVector5_10(k);

//        // 6

//        WVector6_0(k) = getMultiplicationRecurrenceRelation(WVector4_4,WVector4_7,k);               // Added because of the mistake found in the complete transformation matrix
////        WVector6_1(k) = getDivisionRecurrenceRelation(XMatrix.row(3),XMatrix.row(9),WVector6_1,k);
//        WVector6_2(k) = getMultiplicationRecurrenceRelation(WVector4_5,WVector4_15,k);
//        WVector6_3(k) = getMultiplicationRecurrenceRelation(WVector4_6,WVector4_9,k);
//        WVector6_4(k) = getMultiplicationRecurrenceRelation(WVector4_5,WVector4_11,k);
////        WVector6_5(k) = getMultiplicationRecurrenceRelation(WVector4_4,XMatrix.row(16),k);
//        WVector6_5(k) = getMultiplicationRecurrenceRelation(WVector4_4,WVector4_38,k);
//        WVector6_6(k) = getMultiplicationRecurrenceRelation(WVector4_2,(WVector6_2+WVector6_0),k); // Changed because of the mistake found in the complete transformation matrix
//        WVector6_7(k) = getMultiplicationRecurrenceRelation(WVector4_33,WVector6_3,k);
//        WVector6_8(k) = getMultiplicationRecurrenceRelation(WVector4_34,(WVector6_4-WVector6_5),k);

//        UMatrix(6,k) = -standardGravitationalParameter*WVector6_1(k)+WVector6_6(k)-WVector6_7(k)-WVector6_8(k);

        // 7

//        UMatrix(7,k) =-thrust/(tudat::physical_constants::STANDARD_EARTH_GRAVITATIONAL_ACCELERATION*specificImpulse);
        UMatrix(7,k) = 0;

//        // 8

//        WVector8_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(1),XMatrix.row(4),k);
//        WVector8_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(2),XMatrix.row(5),k);
//        WVector8_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(3),XMatrix.row(6),k);

//        UMatrix(8,k) = 2.0*(WVector8_1(k)+WVector8_2(k)+WVector8_3(k));

//        // 9

//        W9IntermediateVector(k) = getMultiplicationRecurrenceRelation(XMatrix.row(9),UMatrix.row(8),k);

//        WVector9(k) = getDivisionRecurrenceRelation(W9IntermediateVector,XMatrix.row(8),WVector9,k);

//        UMatrix(9,k) = 1.5*WVector9(k);

//        // 10

//        UMatrix(10,k) = 0;


        // 11
        WVector11_0(k) = getMultiplicationRecurrenceRelation(XMatrix.row(15),WVector4_38,k);
        WVector11_1(k) = getDivisionRecurrenceRelation(WVector11_0,XMatrix.row(16),WVector11_1,k);

        double Wmult_1 = 0;                   // Setting the output to zero initially (just in case)
        double Wmult_2 = 0;                   // Setting the output to zero initially (just in case)

        for (int j=0; j < k+1; j++){                    // It goes till the order (till k), and stops as soon as j becomes k+1

        if (j == 0){
             Wmult_1 += XMatrix(15,j)*WVector4_38(k-j);          // Wmult += ... means Wmult = Wmult + ...
      }
        else{
            Wmult_1 += XMatrix(15,j)*WVector4_38(k-j);          // Wmult += ... means Wmult = Wmult + ...
            Wmult_2 += XMatrix(16,j)*WVectorTest(k-j);          // Wmult += ... means Wmult = Wmult + ...

        }
        };

        WVectorTest(k) = (Wmult_1-Wmult_2)/XMatrix(16,0);



        WVector11_2(k) = getMultiplicationRecurrenceRelation(WVector11_1,WVector4_9,k);
        WVector11_3(k) = getDivisionRecurrenceRelation(WVector11_2,WVector4_6,WVector11_3,k);

        UMatrix(11,k) = WVector11_3(k);


        // 12
        WVector12_1(k) = getMultiplicationRecurrenceRelation(WVector11_1,WVector4_5,k);

        UMatrix(12,k) = WVector12_1(k);

        // 13
        WVector13_0(k) = getMultiplicationRecurrenceRelation(WVector4_28,WVector4_27,k);
        WVector13_1(k) = getMultiplicationRecurrenceRelation(WVector4_7,WVector4_5,k);
        WVector13_2(k) = rotationalVelocity*rotationalVelocity*getMultiplicationRecurrenceRelation(XMatrix.row(16),WVector4_6,k);
        WVector13_3(k) = getMultiplicationRecurrenceRelation(WVector13_2,WVector4_4,k);
        WVector13_4(k) = Thrust_*getMultiplicationRecurrenceRelation(WVector13_0,WVector4_37,k);
        WVector13_5(k) = getMultiplicationRecurrenceRelation(WVector13_3,WVector4_9,k)+WVector13_4(k);
        WVector13_6(k) = getDivisionRecurrenceRelation(WVector13_5,XMatrix.row(15),WVector13_6,k);
        WVector13_7(k) = -2.0*rotationalVelocity*getMultiplicationRecurrenceRelation(WVector4_6,WVector13_1,k)+WVector13_6(k);
        WVector13_8(k) = getDivisionRecurrenceRelation(WVector13_7,WVector4_38,WVector13_8,k);
        WVector13_9(k) = 2.0*rotationalVelocity*WVector4_4(k)+getMultiplicationRecurrenceRelation(WVector11_3,WVector4_4,k)+WVector13_8(k);

        UMatrix(13,k) = WVector13_9(k);

        // 14
        WVector14_0(k) = Thrust_*getMultiplicationRecurrenceRelation(WVector4_29,WVector4_37,k);
        WVector14_1(k) = getMultiplicationRecurrenceRelation(WVector4_6,WVector4_38,k)+getMultiplicationRecurrenceRelation(WVector13_1,WVector4_4,k);
        WVector14_2(k) = -standardGravitationalParameter*WVector4_38(k);
        WVector14_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(16),XMatrix.row(16),k);
        WVector14_4(k) = getDivisionRecurrenceRelation(WVector14_2,WVector14_3,WVector14_4,k);
        WVector14_5(k) = WVector14_4(k)+getMultiplicationRecurrenceRelation(WVector13_2,WVector14_1,k)+WVector14_0(k);
        WVector14_6(k) = getDivisionRecurrenceRelation(WVector14_5,XMatrix.row(15),WVector14_6,k);
        WVector14_7(k) = 2.0*rotationalVelocity*getMultiplicationRecurrenceRelation(WVector4_6,WVector4_9,k)+WVector11_1(k)+WVector14_6(k);

        UMatrix(14,k) = WVector14_7(k);

        // 15
        WVector15_0(k) = Thrust_*getMultiplicationRecurrenceRelation(WVector4_26,WVector4_27,k)-XMatrix(27,k);
        WVector15_1(k) = getDivisionRecurrenceRelation(WVector15_0,XMatrix.row(7),WVector15_1,k);
//        WVector15_1(k) = getMultiplicationRecurrenceRelation(WVector15_0,WVector4_37,k);
        WVector15_2(k) = -standardGravitationalParameter*WVector4_7(k);
        WVector15_3(k) = getDivisionRecurrenceRelation(WVector15_2,WVector14_3,WVector15_3,k);
        WVector15_4(k) = getMultiplicationRecurrenceRelation(WVector4_38,WVector4_5,k);
        WVector15_5(k) = getMultiplicationRecurrenceRelation(WVector4_7,WVector4_6,k)-getMultiplicationRecurrenceRelation(WVector15_4,WVector4_4,k);
        WVector15_6(k) = getMultiplicationRecurrenceRelation(WVector13_2,WVector15_5,k)+WVector15_1(k)+WVector15_3(k);

        UMatrix(15,k) = WVector15_6(k);

        // 16
        WVector16_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(15),WVector4_7,k);

        UMatrix(16,k) = WVector16_1(k);

        // 31

        UMatrix(31,k) = UMatrix(16,k);



        // 30

        WVector30_1(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_1,9,k);
        WVector30_2(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_2,8,k);
        WVector30_3(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_3,7,k);
        WVector30_4(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_4,6,k);
        WVector30_5(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_5,5,k);
        WVector30_6(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_6,4,k);
        WVector30_7(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_7,3,k);
        WVector30_8(k) = getPowerRecurrenceRelation(XMatrix.row(31),WVector30_8,2,k);

        // Deal with wrong vector orientation...
        Eigen::VectorXd x31Transpose = XMatrix.row(31);
//        std::cout<<"x31Transpose = "<<x31Transpose<<std::endl;
//        x31.row(0)=XMatrix.row(31);
//        Eigen::MatrixXd x31Transpose = x31.transpose();

        WVector30_9(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),(10.0*densityPolyCoefficients(10)*WVector30_1+9.0*densityPolyCoefficients(9)*WVector30_2+8.0*densityPolyCoefficients(8)*WVector30_3+
                                          7.0*densityPolyCoefficients(7)*WVector30_4+6.0*densityPolyCoefficients(6)*WVector30_5+5.0*densityPolyCoefficients(5)*WVector30_6+
                                          4.0*densityPolyCoefficients(4)*WVector30_7+3.0*densityPolyCoefficients(3)*WVector30_8+2.0*densityPolyCoefficients(2)*x31Transpose+densityPolyCoefficient_1),k);

        WVector30_19(k) = (10.0*densityPolyCoefficients(10)*WVector30_1(k)+9.0*densityPolyCoefficients(9)*WVector30_2(k)+8.0*densityPolyCoefficients(8)*WVector30_3(k)+
                7.0*densityPolyCoefficients(7)*WVector30_4(k)+6.0*densityPolyCoefficients(6)*WVector30_5(k)+5.0*densityPolyCoefficients(5)*WVector30_6(k)+
                4.0*densityPolyCoefficients(4)*WVector30_7(k)+3.0*densityPolyCoefficients(3)*WVector30_8(k)+2.0*densityPolyCoefficients(2)*XMatrix(31,k));


        WVector30_9(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),WVector30_19,k) + densityPolyCoefficients(1)*UMatrix(31,k);

        UMatrix(30,k) = WVector30_9(k);



        // 34
        // First it has to be determined which equation has to be used and which values have to be computed

        if (temperatureAltitudeRanges(0,0)<=XMatrix(31,0) && XMatrix(31,0) < temperatureAltitudeRanges(0,1)){

            UMatrix(34,k) = temperaturePolyCoefficients(0,1)*UMatrix(31,k);
        }

        else if (temperatureAltitudeRanges(1,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(1,1)){

            WVector34_12(k) = (2.0*temperaturePolyCoefficients(1,2)*XMatrix(31,k));

            WVector34_2(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),WVector34_12,k)+temperaturePolyCoefficients(1,1)*UMatrix(31,k);
        // change to Eigen::VectorXd

            UMatrix(34,k) = WVector34_2(k);
        }

        else if (temperatureAltitudeRanges(2,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(2,1)){

//            WVector34_3(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),(6.0*temperaturePolyCoefficients(2,6)*WVector30_5+5.0*temperaturePolyCoefficients(2,5)*WVector30_6+4.0*temperaturePolyCoefficients(2,4)*WVector30_7+
//                                                                                  3.0*temperaturePolyCoefficients(2,3)*WVector30_8+2.0*temperaturePolyCoefficients(2,2)*x31Transpose+temperaturePolyCoefficient_2_1),k);

            WVector34_13(k) = (6.0*temperaturePolyCoefficients(2,6)*WVector30_5(k)+5.0*temperaturePolyCoefficients(2,5)*WVector30_6(k)+4.0*temperaturePolyCoefficients(2,4)*WVector30_7(k)+
                               3.0*temperaturePolyCoefficients(2,3)*WVector30_8(k)+2.0*temperaturePolyCoefficients(2,2)*XMatrix(31,k));

            WVector34_3(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),WVector34_13,k)+temperaturePolyCoefficients(2,1)*UMatrix(31,k);



            UMatrix(34,k) = WVector34_3(k);

            }

        else if (temperatureAltitudeRanges(3,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(3,1)){

//            WVector34_4(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),(8.0*temperaturePolyCoefficients(3,8)*WVector30_3+7.0*temperaturePolyCoefficients(3,7)*WVector30_4+
//                                                                     6.0*temperaturePolyCoefficients(3,6)*WVector30_5+5.0*temperaturePolyCoefficients(3,5)*WVector30_6+4.0*temperaturePolyCoefficients(3,4)*WVector30_7+
//                                                                                  3.0*temperaturePolyCoefficients(3,3)*WVector30_8+2.0*temperaturePolyCoefficients(3,2)*x31Transpose+temperaturePolyCoefficient_3_1),k);

            WVector34_14(k) = (8.0*temperaturePolyCoefficients(3,8)*WVector30_3(k)+7.0*temperaturePolyCoefficients(3,7)*WVector30_4(k)+
                               6.0*temperaturePolyCoefficients(3,6)*WVector30_5(k)+5.0*temperaturePolyCoefficients(3,5)*WVector30_6(k)+4.0*temperaturePolyCoefficients(3,4)*WVector30_7(k)+
                                            3.0*temperaturePolyCoefficients(3,3)*WVector30_8(k)+2.0*temperaturePolyCoefficients(3,2)*XMatrix(31,k));

            WVector34_4(k) = getMultiplicationRecurrenceRelation(UMatrix.row(31),WVector34_14,k)+temperaturePolyCoefficients(3,1)*UMatrix(31,k);

            UMatrix(34,k) = WVector34_4(k);


        }
        else if (temperatureAltitudeRanges(4,0)<=XMatrix(31,0) && XMatrix(31,0)<temperatureAltitudeRanges(4,1)){

            UMatrix(34,k) = 0;


    };



        // 28

        WVector28(k) = getMultiplicationRecurrenceRelation(UMatrix.row(30),XMatrix.row(28),k);

        UMatrix(28,k) = WVector28(k);

        // 33

        WVector33(k) = getDivisionRecurrenceRelation(UMatrix.row(34),XMatrix.row(33),WVector33,k);

        UMatrix(33,k) = WVector33(k)*((adiabeticIndex*specificGasConstant)/2.0);


        // 32

        WVector32_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(33),UMatrix.row(15),k);
        WVector32_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(15),UMatrix.row(33),k);
        WVector32_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(33),XMatrix.row(33),k);
        WVector32_4(k) = getDivisionRecurrenceRelation((WVector32_1-WVector32_2),WVector32_3,WVector32_4,k);

        UMatrix(32,k) = WVector32_4(k);

        // 29
        // Determine which section of the drag coefficient curve needs to be used

        for (int i=0; i < 5+1; i++){

            if (dragCoefficientMachRanges(i,0) <= XMatrix(32,0) && XMatrix(32,0) < dragCoefficientMachRanges(i,1)){

                sectionCD = i;
            }

        };


        UMatrix(29,k) = dragCoefficientPolyCoefficients(sectionCD,1)*UMatrix(32,k);


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

    stateTaylorCoefficients.row(1) = XMatrix.row(16);   // Radius (r)  [km]
    stateTaylorCoefficients.row(2) = XMatrix.row(12);   // Latitude (delta)  [rad]
    stateTaylorCoefficients.row(3) = XMatrix.row(11);   // Longitude (tau)  [rad]
    stateTaylorCoefficients.row(4) = XMatrix.row(15);   // Ground velocity (V_G) [km/s]
    stateTaylorCoefficients.row(5) = XMatrix.row(14);   // Flight-path angle (gamma_G)  [rad]
    stateTaylorCoefficients.row(6) = XMatrix.row(13);   // Azimuth angle (chi_G)    [rad]
    stateTaylorCoefficients.row(7) = XMatrix.row(7);    // Mass [kg] from literature study



    stateTaylorCoefficients.row(0) = XMatrix.row(32);   // Mach number [-] required to determine the time at which the CD section boundary was passsed

    /// Debug ///

//    std::cout<<"WVector4_7 = "<<WVector4_7<<std::endl;


//    std::cout<<"WVector11_1 = "<<WVector11_1<<std::endl;
//    std::cout<<"WVector4_12 = "<<WVector4_12<<std::endl;

//    std::cout<<"WVector13_0 = "<<WVector13_0<<std::endl;
//    std::cout<<"WVector13_1 = "<<WVector13_1<<std::endl;
//    std::cout<<"WVector13_2 = "<<WVector13_2<<std::endl;
//    std::cout<<"WVector13_3 = "<<WVector13_3<<std::endl;
//    std::cout<<"WVector13_4 = "<<WVector13_4<<std::endl;
//    std::cout<<"WVector13_5 = "<<WVector13_5<<std::endl;
//    std::cout<<"WVector13_6 = "<<WVector13_6<<std::endl;
//    std::cout<<"WVector13_7 = "<<WVector13_7<<std::endl;
//    std::cout<<"WVector13_8 = "<<WVector13_8<<std::endl;
//    std::cout<<"WVector13_9 = "<<WVector13_9<<std::endl;

//    std::cout<<"WVector13_10 = "<<WVector13_10<<std::endl;
//    std::cout<<"WVector13_11 = "<<WVector13_11<<std::endl;

//    std::cout<<"WVector4_6 (cos(x12)) = "<<WVector4_6<<std::endl;

//    std::cout<<"WVector4_3 (sin(x10+x11)) = "<<WVector4_3<<std::endl;
//    std::cout<<"WVector4_8 (cos(x10+x11)) = "<<WVector4_8<<std::endl;

//    std::cout<<"WVector4_5 (cos(x13)) = "<<WVector4_5<<std::endl;
//    std::cout<<"WVector4_9 (sin(x13)) = "<<WVector4_9<<std::endl;

//    std::cout<<"WVector4_7 (sin(x14)) = "<<WVector4_7<<std::endl;
//    std::cout<<"WVector4_38 (cos(x14)) = "<<WVector4_38<<std::endl;

//    std::cout<<"WVector13_1 (sin(x14)*cos(x13)) = "<<WVector13_1<<std::endl;
//    std::cout<<"WVector4_7 (sin(x14)) = "<<WVector4_7<<std::endl;
//    std::cout<<"WVector4_5 (cos(x13)) = "<<WVector4_5<<std::endl;
//    std::cout<<"WVector4_4 (sin(x12)) = "<<WVector4_4<<std::endl;

//    std::cout<<"WVector14_0 = "<<WVector14_0<<std::endl;
//    std::cout<<"WVector14_1 = "<<WVector14_1<<std::endl;
//    std::cout<<"WVector14_2 = "<<WVector14_2<<std::endl;
//    std::cout<<"WVector14_3 = "<<WVector14_3<<std::endl;
//    std::cout<<"WVector14_4 = "<<WVector14_4<<std::endl;
//    std::cout<<"WVector14_5 = "<<WVector14_5<<std::endl;
//    std::cout<<"WVector14_6 = "<<WVector14_6<<std::endl;
//    std::cout<<"WVector14_7 = "<<WVector14_7<<std::endl;
//    std::cout<<"WVector14_8 = "<<WVector14_8<<std::endl;
//    std::cout<<"WVector14_9 = "<<WVector14_9<<std::endl;

//    std::cout<<"WVector11_1 = "<<WVector11_1<<std::endl;
//    std::cout<<"WVectorTest = "<<WVectorTest<<std::endl;
//    std::cout<<"WVector11_2 = "<<WVector11_2<<std::endl;
//    std::cout<<"WVector11_3 = "<<WVector11_3<<std::endl;
//    std::cout<<"WVector4_6 = "<<WVector4_6<<std::endl;
//    std::cout<<"WVector4_9 = "<<WVector4_9<<std::endl;


//        std::cout<<"WVector15_0 = "<<WVector15_0<<std::endl;
//        std::cout<<"WVector15_1 = "<<WVector15_1<<std::endl;
//        std::cout<<"WVector15_6 = "<<WVector15_6<<std::endl;
//        std::cout<<"WVector4_25 = "<<WVector4_25<<std::endl;
//        std::cout<<"WVector15_2 = "<<WVector15_2<<std::endl;
//        std::cout<<"WVector15_3 = "<<WVector15_3<<std::endl;
//        std::cout<<"WVector15_4 = "<<WVector15_4<<std::endl;
//        std::cout<<"WVector15_5 = "<<WVector15_5<<std::endl;
//        std::cout<<"WVector15_6 = "<<WVector15_6<<std::endl;



//    std::cout<<"WVector4_6 = "<<WVector4_6<<std::endl;
//    std::cout<<"WVector4_9 = "<<WVector4_9<<std::endl;
//    std::cout<<"WVector4_38 = "<<WVector4_38<<std::endl;

//    std::cout<<"XMatrix.row(4) = "<<XMatrix.row(4)<<std::endl;
//    std::cout<<"XMatrix.row(5) = "<<XMatrix.row(5)<<std::endl;

//std::cout<<"XMatrix.row(11) = "<<XMatrix.row(11)<<std::endl;
//std::cout<<"XMatrix.row(12) = "<<XMatrix.row(12)<<std::endl;
//std::cout<<"XMatrix.row(13) = "<<XMatrix.row(13)<<std::endl;
//std::cout<<"XMatrix.row(14) = "<<XMatrix.row(14)<<std::endl;
//std::cout<<"XMatrix.row(15) = "<<XMatrix.row(15)<<std::endl;
//std::cout<<"XMatrix.row(16) = "<<XMatrix.row(16)<<std::endl;

//std::cout<<"UMatrix.row(11) = "<<UMatrix.row(11)<<std::endl;
//std::cout<<"UMatrix.row(12) = "<<UMatrix.row(12)<<std::endl;
//std::cout<<"UMatrix.row(13) = "<<UMatrix.row(13)<<std::endl;
//std::cout<<"UMatrix.row(14) = "<<UMatrix.row(14)<<std::endl;
//std::cout<<"UMatrix.row(15) = "<<UMatrix.row(15)<<std::endl;
//std::cout<<"UMatrix.row(16) = "<<UMatrix.row(16)<<std::endl;

//  std::cout<<"XMatrix.row(8) = "<<XMatrix.row(8)<<std::endl;
//  std::cout<<"XMatrix.row(9) = "<<XMatrix.row(9)<<std::endl;

//    std::cout<<"WVector6_0 = "<<WVector6_0<<std::endl;
//    std::cout<<"WVector6_1 = "<<WVector6_1<<std::endl;
//    std::cout<<"WVector6_2 = "<<WVector6_2<<std::endl;
//    std::cout<<"WVector6_3 = "<<WVector6_3<<std::endl;
//    std::cout<<"WVector6_4 = "<<WVector6_4<<std::endl;
//    std::cout<<"WVector6_5 = "<<WVector6_5<<std::endl;
//    std::cout<<"WVector6_6 = "<<WVector6_6<<std::endl;
//    std::cout<<"WVector6_7 = "<<WVector6_7<<std::endl;
//    std::cout<<"WVector6_8 = "<<WVector6_8<<std::endl;

//std::cout<<"WVector32_1 = "<<WVector32_1<<std::endl;
//std::cout<<"WVector32_2 = "<<WVector32_2<<std::endl;
//std::cout<<"WVector32_3 = "<<WVector32_3<<std::endl;
//std::cout<<"WVector32_4 = "<<WVector32_4<<std::endl;

//    std::cout<<"WVector34_14 = "<<WVector34_14<<std::endl;
//    std::cout<<"WVector34_4 = "<<WVector34_4<<std::endl;

//    std::cout<<"WVector4_26 = "<<WVector4_26<<std::endl;
//    std::cout<<"WVector4_28 = "<<WVector4_28<<std::endl;

//    std::cout<<"WVector27_1 = "<<WVector27_1<<std::endl;
//    std::cout<<"WVector27_2 = "<<WVector27_2<<std::endl;
//    std::cout<<"WVector27_3 = "<<WVector27_3<<std::endl;
//    std::cout<<"WVector27_4 = "<<WVector27_4<<std::endl;
//    std::cout<<"WVector27_5 = "<<WVector27_5<<std::endl;
//    std::cout<<"WVector27_6 = "<<WVector27_6<<std::endl;



//std::cout<<"XMatrix.row(33) = "<<XMatrix.row(33)<<std::endl;
//std::cout<<"XMatrix.row(34) = "<<XMatrix.row(34)<<std::endl;
//std::cout<<"XMatrix.row(15) = "<<XMatrix.row(15)<<std::endl;
//std::cout<<"XMatrix.row(31) = "<<XMatrix.row(31)<<std::endl;
//std::cout<<"XMatrix.row(16) = "<<XMatrix.row(16)<<std::endl;
//    std::cout<<"XMatrix.row(27) = "<<XMatrix.row(27)<<std::endl;
//    std::cout<<"XMatrix.row(28) = "<<XMatrix.row(28)<<std::endl;


//std::cout<<"XMatrix.row(32) = "<<XMatrix.row(32)<<std::endl;


//    std::cout<<"XMatrix = "<<XMatrix<<std::endl;
//    std::cout<<"UMatrix = "<<UMatrix<<std::endl;






    return stateTaylorCoefficients;

} // end of function




