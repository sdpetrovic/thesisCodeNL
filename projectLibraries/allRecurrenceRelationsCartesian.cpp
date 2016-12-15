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



#include "allRecurrenceRelationsCartesian.h"
#include "basicRecurrenceRelations.h"



Eigen::MatrixXd getCartesianTaylorCoefficients(const double adiabeticIndex_, const double specificGasConstant_, const double standardGravitationalParameter_, const double rotationalVelocity_, const double primeMeridianAngle_,
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

//    Eigen::VectorXd WVector4_0 = Eigen::VectorXd::Zero(maxOrder);     // W4,0   // Added because of the mistake found in the recurrence relation of W4,2
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
    Eigen::VectorXd WVector4_39 = Eigen::VectorXd::Zero(maxOrder);     // W4,39
    Eigen::VectorXd WVector4_40 = Eigen::VectorXd::Zero(maxOrder);     // W4,40

    Eigen::VectorXd WVector4_41 = Eigen::VectorXd::Zero(maxOrder);     // W4,41
    Eigen::VectorXd WVector4_42 = Eigen::VectorXd::Zero(maxOrder);     // W4,42
    Eigen::VectorXd WVector4_43 = Eigen::VectorXd::Zero(maxOrder);     // W4,43
    Eigen::VectorXd WVector4_44 = Eigen::VectorXd::Zero(maxOrder);     // W4,44
    Eigen::VectorXd WVector4_45 = Eigen::VectorXd::Zero(maxOrder);     // W4,45
    Eigen::VectorXd WVector4_46 = Eigen::VectorXd::Zero(maxOrder);     // W4,46
    Eigen::VectorXd WVector4_47 = Eigen::VectorXd::Zero(maxOrder);     // W4,47
    Eigen::VectorXd WVector4_48 = Eigen::VectorXd::Zero(maxOrder);     // W4,48
    Eigen::VectorXd WVector4_49 = Eigen::VectorXd::Zero(maxOrder);     // W4,49
    Eigen::VectorXd WVector4_50 = Eigen::VectorXd::Zero(maxOrder);     // W4,50

    Eigen::VectorXd WVector4_51 = Eigen::VectorXd::Zero(maxOrder);     // W4,51
    Eigen::VectorXd WVector4_52 = Eigen::VectorXd::Zero(maxOrder);     // W4,52


    Eigen::VectorXd WVector5_1 = Eigen::VectorXd::Zero(maxOrder);     // W5,1
    Eigen::VectorXd WVector5_2 = Eigen::VectorXd::Zero(maxOrder);     // W5,2
    Eigen::VectorXd WVector5_3 = Eigen::VectorXd::Zero(maxOrder);     // W5,3
    Eigen::VectorXd WVector5_4 = Eigen::VectorXd::Zero(maxOrder);     // W5,4
    Eigen::VectorXd WVector5_5 = Eigen::VectorXd::Zero(maxOrder);     // W5,5
    Eigen::VectorXd WVector5_6 = Eigen::VectorXd::Zero(maxOrder);     // W5,6
//    Eigen::VectorXd WVector5_7 = Eigen::VectorXd::Zero(maxOrder);     // W5,7
//    Eigen::VectorXd WVector5_8 = Eigen::VectorXd::Zero(maxOrder);     // W5,8
//    Eigen::VectorXd WVector5_9 = Eigen::VectorXd::Zero(maxOrder);     // W5,9
//    Eigen::VectorXd WVector5_10 = Eigen::VectorXd::Zero(maxOrder);     // W5,10


//    Eigen::VectorXd WVector6_0 = Eigen::VectorXd::Zero(maxOrder);     // W6,0  // Added because of the mistake found in the complete transformation matrix
    Eigen::VectorXd WVector6_1 = Eigen::VectorXd::Zero(maxOrder);     // W6,1
    Eigen::VectorXd WVector6_2 = Eigen::VectorXd::Zero(maxOrder);     // W6,2
    Eigen::VectorXd WVector6_3 = Eigen::VectorXd::Zero(maxOrder);     // W6,3
    Eigen::VectorXd WVector6_4 = Eigen::VectorXd::Zero(maxOrder);     // W6,4
    Eigen::VectorXd WVector6_5 = Eigen::VectorXd::Zero(maxOrder);     // W6,5
    Eigen::VectorXd WVector6_6 = Eigen::VectorXd::Zero(maxOrder);     // W6,6
    Eigen::VectorXd WVector6_7 = Eigen::VectorXd::Zero(maxOrder);     // W6,7
//    Eigen::VectorXd WVector6_8 = Eigen::VectorXd::Zero(maxOrder);     // W6,8


    Eigen::VectorXd WVector8_1 = Eigen::VectorXd::Zero(maxOrder);     // W8,1
    Eigen::VectorXd WVector8_2 = Eigen::VectorXd::Zero(maxOrder);     // W8,2
    Eigen::VectorXd WVector8_3 = Eigen::VectorXd::Zero(maxOrder);     // W8,3


    Eigen::VectorXd WVector9 = Eigen::VectorXd::Zero(maxOrder);     // W9


    Eigen::VectorXd WVector27_1 = Eigen::VectorXd::Zero(maxOrder);     // W27,1
    Eigen::VectorXd WVector27_2 = Eigen::VectorXd::Zero(maxOrder);     // W27,2
    Eigen::VectorXd WVector27_3 = Eigen::VectorXd::Zero(maxOrder);     // W27,3
    Eigen::VectorXd WVector27_4 = Eigen::VectorXd::Zero(maxOrder);     // W27,4
    Eigen::VectorXd WVector27_5 = Eigen::VectorXd::Zero(maxOrder);     // W27,5
    Eigen::VectorXd WVector27_6 = Eigen::VectorXd::Zero(maxOrder);     // W27,6
    Eigen::VectorXd WVector27_7 = Eigen::VectorXd::Zero(maxOrder);     // W27,7
    Eigen::VectorXd WVector27_8 = Eigen::VectorXd::Zero(maxOrder);     // W27,8
    Eigen::VectorXd WVector27_9 = Eigen::VectorXd::Zero(maxOrder);     // W27,9
    Eigen::VectorXd WVector27_10 = Eigen::VectorXd::Zero(maxOrder);     // W27,10

    Eigen::VectorXd WVector27_11 = Eigen::VectorXd::Zero(maxOrder);     // W27,11
    Eigen::VectorXd WVector27_12 = Eigen::VectorXd::Zero(maxOrder);     // W27,12
    Eigen::VectorXd WVector27_13 = Eigen::VectorXd::Zero(maxOrder);     // W27,13
    Eigen::VectorXd WVector27_14 = Eigen::VectorXd::Zero(maxOrder);     // W27,14
    Eigen::VectorXd WVector27_15 = Eigen::VectorXd::Zero(maxOrder);     // W27,15
    Eigen::VectorXd WVector27_16 = Eigen::VectorXd::Zero(maxOrder);     // W27,16
    Eigen::VectorXd WVector27_17 = Eigen::VectorXd::Zero(maxOrder);     // W27,17
    Eigen::VectorXd WVector27_18 = Eigen::VectorXd::Zero(maxOrder);     // W27,18
    Eigen::VectorXd WVector27_19 = Eigen::VectorXd::Zero(maxOrder);     // W27,19






    // And fill them

//    WVector4_0(0) = initialFunctionsMatrix(4,0);     // W4,0   // Added because of the mistake found in the recurrence relation of W4,2
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
    WVector4_39(0) = initialFunctionsMatrix(4,39);     // W4,39
    WVector4_40(0) = initialFunctionsMatrix(4,40);     // W4,40

    WVector4_41(0) = initialFunctionsMatrix(4,41);     // W4,41
    WVector4_42(0) = initialFunctionsMatrix(4,42);     // W4,42
    WVector4_43(0) = initialFunctionsMatrix(4,43);     // W4,43
    WVector4_44(0) = initialFunctionsMatrix(4,44);     // W4,44
    WVector4_45(0) = initialFunctionsMatrix(4,45);     // W4,45
    WVector4_46(0) = initialFunctionsMatrix(4,46);     // W4,46
    WVector4_47(0) = initialFunctionsMatrix(4,47);     // W4,47
    WVector4_48(0) = initialFunctionsMatrix(4,48);     // W4,48
    WVector4_49(0) = initialFunctionsMatrix(4,49);     // W4,49
    WVector4_50(0) = initialFunctionsMatrix(4,50);     // W4,50

    WVector4_51(0) = initialFunctionsMatrix(4,51);     // W4,51
    WVector4_52(0) = initialFunctionsMatrix(4,52);     // W4,52

    WVector5_1(0) = initialFunctionsMatrix(5,1);     // W5,1
    WVector5_2(0) = initialFunctionsMatrix(5,2);     // W5,2
    WVector5_3(0) = initialFunctionsMatrix(5,3);     // W5,3
    WVector5_4(0) = initialFunctionsMatrix(5,4);     // W5,4
    WVector5_5(0) = initialFunctionsMatrix(5,5);     // W5,5
    WVector5_6(0) = initialFunctionsMatrix(5,6);     // W5,6
//    WVector5_7(0) = initialFunctionsMatrix(5,7);     // W5,7
//    WVector5_8(0) = initialFunctionsMatrix(5,8);     // W5,8
//    WVector5_9(0) = initialFunctionsMatrix(5,9);     // W5,9
//    WVector5_10(0) = initialFunctionsMatrix(5,10);     // W5,10


//    WVector6_0(0) = initialFunctionsMatrix(6,0);     // W6,0  // Added because of the mistake found in the complete transformation matrix
    WVector6_1(0) = initialFunctionsMatrix(6,1);     // W6,1
    WVector6_2(0) = initialFunctionsMatrix(6,2);     // W6,2
    WVector6_3(0) = initialFunctionsMatrix(6,3);     // W6,3
    WVector6_4(0) = initialFunctionsMatrix(6,4);     // W6,4
    WVector6_5(0) = initialFunctionsMatrix(6,5);     // W6,5
    WVector6_6(0) = initialFunctionsMatrix(6,6);     // W6,6
    WVector6_7(0) = initialFunctionsMatrix(6,7);     // W6,7
//    WVector6_8(0) = initialFunctionsMatrix(6,8);     // W6,8


    WVector8_1(0) = initialFunctionsMatrix(8,1);     // W8,1
    WVector8_2(0) = initialFunctionsMatrix(8,2);     // W8,2
    WVector8_3(0) = initialFunctionsMatrix(8,3);     // W8,3


    WVector9(0) = initialFunctionsMatrix(9,1);     // W9




    WVector27_1(0) = initialFunctionsMatrix(27,1);     // W27,1
    WVector27_2(0) = initialFunctionsMatrix(27,2);     // W27,2
    WVector27_3(0) = initialFunctionsMatrix(27,3);     // W27,3
    WVector27_4(0) = initialFunctionsMatrix(27,4);     // W27,4
    WVector27_5(0) = initialFunctionsMatrix(27,5);     // W27,5
    WVector27_6(0) = initialFunctionsMatrix(27,6);     // W27,6
    WVector27_7(0) = initialFunctionsMatrix(27,7);     // W27,7
    WVector27_8(0) = initialFunctionsMatrix(27,8);     // W27,8
    WVector27_9(0) = initialFunctionsMatrix(27,9);     // W27,9
    WVector27_10(0) = initialFunctionsMatrix(27,10);     // W27,10

    WVector27_11(0) = initialFunctionsMatrix(27,11);     // W27,11
    WVector27_12(0) = initialFunctionsMatrix(27,12);     // W27,12
    WVector27_13(0) = initialFunctionsMatrix(27,13);     // W27,13
    WVector27_14(0) = initialFunctionsMatrix(27,14);     // W27,14
    WVector27_15(0) = initialFunctionsMatrix(27,15);     // W27,15
    WVector27_16(0) = initialFunctionsMatrix(27,16);     // W27,16
    WVector27_17(0) = initialFunctionsMatrix(27,17);     // W27,17
    WVector27_18(0) = initialFunctionsMatrix(27,18);     // W27,18
    WVector27_19(0) = initialFunctionsMatrix(27,19);     // W27,19




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
    Eigen::VectorXd temperaturePolyCoefficient_1_2 = Eigen::VectorXd::Zero(maxOrder);
    Eigen::VectorXd temperaturePolyCoefficient_1_3 = Eigen::VectorXd::Zero(maxOrder);
    Eigen::VectorXd temperaturePolyCoefficient_1_4 = Eigen::VectorXd::Zero(maxOrder);
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
            temperaturePolyCoefficient_1_2(i) = temperaturePolyCoefficients(1,0);
            temperaturePolyCoefficient_1_3(i) = temperaturePolyCoefficients(2,0);
            temperaturePolyCoefficient_1_4(i) = temperaturePolyCoefficients(3,0);
            ThrustVector(i) = Thrust_;
            if (i==0){
            thrustAzimuthVector(i) = thrustAzimuthMatrix(0,2);
            thrustElevationVector(i) = thrustElevationMatrix(0,2);

//            onesVector(i) = 1;
};

        }


//    std::cout<<"So far so good 3"<<std::endl;

        // 1
        UMatrix(1,k) = XMatrix(4,k);

        // 2
        UMatrix(2,k) = XMatrix(5,k);

        // 3
        UMatrix(3,k) = XMatrix(6,k);

        // 4

        WVector4_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(1),XMatrix.row(1),k)+getMultiplicationRecurrenceRelation(XMatrix.row(2),XMatrix.row(2),k);
        WVector4_2(k) = WVector4_1(k)+getMultiplicationRecurrenceRelation(XMatrix.row(3),XMatrix.row(3),k);

//        if (WVector4_1 == Eigen::VectorXd::Zero(k+1)){
//            WVector4_3(k) = fabs(XMatrix(3,k));
//        }
//        else if (XMatrix.row(1) == Eigen::VectorXd::Zero(k+1) && XMatrix.row(3) == Eigen::VectorXd::Zero(k+1) ){
//            WVector4_3(k) = fabs(XMatrix(2,k));
//        }
//        else if (XMatrix.row(2) == Eigen::VectorXd::Zero(k+1) && XMatrix.row(3) == Eigen::VectorXd::Zero(k+1) ){
//            WVector4_3(k) = fabs(XMatrix(1,k));
////            std::cout<<"It doesn't actually go here does it?"<<std::endl;
//        }
//        else {
        WVector4_3(k) = getPowerRecurrenceRelation(WVector4_2,WVector4_3,0.5,k); // Radius
//        }

//        if (XMatrix.row(1) == Eigen::VectorXd::Zero(k+1)){
//            WVector4_4(k) = fabs(XMatrix(2,k));
//        }
//        else if (XMatrix.row(2) == Eigen::VectorXd::Zero(k+1)){
//            WVector4_4(k) = fabs(XMatrix(1,k));
//        }
//        else{
        WVector4_4(k) = getPowerRecurrenceRelation(WVector4_1,WVector4_4,0.5,k);
//        }
//        if (k == 1 && XMatrix(2,1)-WVector4_4(1)*WVector4_5(0) <= 1e-18){
//            WVector4_5(k) = 0.0;
//        }
//        else{
        WVector4_5(k) = getDivisionRecurrenceRelation(XMatrix.row(2),WVector4_4,WVector4_5,k);  // sin(lambda)
//        }
        WVector4_6(k) = getDivisionRecurrenceRelation(XMatrix.row(1),WVector4_4,WVector4_6,k);  // cos(lambda)
        WVector4_7(k) = getDivisionRecurrenceRelation(XMatrix.row(3),WVector4_3,WVector4_7,k);  // sin(delta)
        WVector4_8(k) = getDivisionRecurrenceRelation(WVector4_4,WVector4_3,WVector4_8,k);      // cos(delta)
        WVector4_9(k) = (XMatrix(4,k)+rotationalVelocity*XMatrix(2,k));
        WVector4_10(k) = (XMatrix(5,k)-rotationalVelocity*XMatrix(1,k));

        WVector4_11(k) = getMultiplicationRecurrenceRelation(WVector4_9,WVector4_9,k)+getMultiplicationRecurrenceRelation(WVector4_10,WVector4_10,k)+getMultiplicationRecurrenceRelation(XMatrix.row(6),XMatrix.row(6),k);

        /// Debug ///
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_9,WVector4_9,k) = "<<getMultiplicationRecurrenceRelation(WVector4_9,WVector4_9,k)<<std::endl;
//        Eigen::Vector3d test1;
//        Eigen::Vector3d test2;
//        test1 << 0, 1, 2;
//        test2 << 0, 1, 4;

//        std::cout<<"getMultiplicationRecurrenceRelation(test1,test2,2) = "<<getMultiplicationRecurrenceRelation(test1,test2,2)<<std::endl;
//        std::cout<<"getPowerRecurrenceRelation(test1,test2,2.0,1.0) = "<<getPowerRecurrenceRelation(test1,test2,2.0,2)<<std::endl;
        /// Debug ///

//        if (WVector4_9 == Eigen::VectorXd::Zero(k+1) && WVector4_10 == Eigen::VectorXd::Zero(k+1)){
//            WVector4_12(k) = fabs(XMatrix(6,k));
//        }
//        else if (WVector4_9 == Eigen::VectorXd::Zero(k+1) && XMatrix.row(6) == Eigen::VectorXd::Zero(k+1) ){
//            WVector4_12(k) = fabs(WVector4_10(k));
//        }
//        else if (WVector4_10 == Eigen::VectorXd::Zero(k+1) && XMatrix.row(6) == Eigen::VectorXd::Zero(k+1) ){
//            WVector4_12(k) = fabs(WVector4_9(k));
////            std::cout<<"It doesn't actually go here does it?"<<std::endl;
//        }
//        else {
         WVector4_12(k) = getPowerRecurrenceRelation(WVector4_11,WVector4_12,0.5,k);    // Ground velocity
//        }



         WVector4_13(k) = getDivisionRecurrenceRelation(WVector4_9,WVector4_12,WVector4_13,k);  //X1
         WVector4_14(k) = getDivisionRecurrenceRelation(WVector4_10,WVector4_12,WVector4_14,k);  //X2
         WVector4_15(k) = getDivisionRecurrenceRelation(XMatrix.row(6),WVector4_12,WVector4_15,k);  //X3

//         WVector4_16(k) = getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(3),k)-getMultiplicationRecurrenceRelation(XMatrix.row(6),XMatrix.row(2),k); // Bergsma
//         WVector4_17(k) = getMultiplicationRecurrenceRelation(XMatrix.row(6),XMatrix.row(1),k)-getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(3),k);
//         WVector4_18(k) = getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(2),k)-getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(1),k);

         WVector4_16(k) = getMultiplicationRecurrenceRelation(XMatrix.row(6),XMatrix.row(2),k)-getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(3),k);
         WVector4_17(k) = getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(3),k)-getMultiplicationRecurrenceRelation(XMatrix.row(6),XMatrix.row(1),k);
         WVector4_18(k) = getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(1),k)-getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(2),k);

//         std::cout<<"getMultiplicationRecurrenceRelation(XMatrix.row(6),XMatrix.row(2),k) = "<<getMultiplicationRecurrenceRelation(XMatrix.row(6),XMatrix.row(2),k)<<std::endl;
//         std::cout<<"-getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(3),k) = "<<-getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(3),k)<<std::endl;
//         std::cout<<"---"<<std::endl;
//         std::cout<<"getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(3),k) = "<<getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(3),k)<<std::endl;
//         std::cout<<"-getMultiplicationRecurrenceRelation(XMatrix.row(6),XMatrix.row(1),k) = "<<-getMultiplicationRecurrenceRelation(XMatrix.row(6),XMatrix.row(1),k)<<std::endl;
//         std::cout<<"---"<<std::endl;
//         std::cout<<"getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(1),k) = "<<getMultiplicationRecurrenceRelation(WVector4_10,XMatrix.row(1),k)<<std::endl;
//         std::cout<<"-getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(2),k) = "<<-getMultiplicationRecurrenceRelation(WVector4_9,XMatrix.row(2),k)<<std::endl;
//         std::cout<<" "<<std::endl;



         WVector4_19(k) = getMultiplicationRecurrenceRelation(WVector4_16,WVector4_16,k)+getMultiplicationRecurrenceRelation(WVector4_17,WVector4_17,k)+getMultiplicationRecurrenceRelation(WVector4_18,WVector4_18,k);
         WVector4_20(k) = getPowerRecurrenceRelation(WVector4_19,WVector4_20,0.5,k);

         WVector4_21(k) = getDivisionRecurrenceRelation(WVector4_16,WVector4_20,WVector4_21,k); // Y1
         WVector4_22(k) = getDivisionRecurrenceRelation(WVector4_17,WVector4_20,WVector4_21,k); // Y2
         WVector4_23(k) = getDivisionRecurrenceRelation(WVector4_18,WVector4_20,WVector4_21,k); // Y3


         WVector4_24(k) = getMultiplicationRecurrenceRelation(WVector4_14,WVector4_23,k)-getMultiplicationRecurrenceRelation(WVector4_15,WVector4_22,k); // Z1
         WVector4_25(k) = getMultiplicationRecurrenceRelation(WVector4_15,WVector4_21,k)-getMultiplicationRecurrenceRelation(WVector4_13,WVector4_23,k); // Z2
         WVector4_40(k) = getMultiplicationRecurrenceRelation(WVector4_13,WVector4_22,k)-getMultiplicationRecurrenceRelation(WVector4_14,WVector4_21,k); // Z3



//        WVector4_13(k) = -getMultiplicationRecurrenceRelation(WVector4_6,WVector4_7,k);
//        WVector4_14(k) = -getMultiplicationRecurrenceRelation(WVector4_7,WVector4_5,k);
//        WVector4_15(k) = -getMultiplicationRecurrenceRelation(WVector4_8,WVector4_6,k);
//        WVector4_16(k) = -getMultiplicationRecurrenceRelation(WVector4_8,WVector4_5,k);
//        WVector4_17(k) = getMultiplicationRecurrenceRelation(XMatrix.row(6),WVector4_8,k)+getMultiplicationRecurrenceRelation(WVector4_9,WVector4_13,k)+getMultiplicationRecurrenceRelation(WVector4_10,WVector4_14,k); // Vx_v
////        WVector4_17(k) = 0.0;

//        /// Debug ///

//        if (k == 1){
//        std::cout<<"getMultiplicationRecurrenceRelation(XMatrix.row(6),WVector4_8,k) = "<<getMultiplicationRecurrenceRelation(XMatrix.row(6),WVector4_8,k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_9,WVector4_13,k) = "<<getMultiplicationRecurrenceRelation(WVector4_9,WVector4_13,k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_10,WVector4_14,k) = "<<getMultiplicationRecurrenceRelation(WVector4_10,WVector4_14,k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(XMatrix.row(6),WVector4_8,k)+getMultiplicationRecurrenceRelation(WVector4_9,WVector4_13,k) = "<<getMultiplicationRecurrenceRelation(XMatrix.row(6),WVector4_8,k)+getMultiplicationRecurrenceRelation(WVector4_9,WVector4_13,k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(XMatrix.row(6),WVector4_8,k)/getMultiplicationRecurrenceRelation(WVector4_9,WVector4_13,k) = "<<getMultiplicationRecurrenceRelation(XMatrix.row(6),WVector4_8,k)/getMultiplicationRecurrenceRelation(WVector4_9,WVector4_13,k)<<std::endl;

//        }
//        /// Debug ///

//        WVector4_18(k) = getMultiplicationRecurrenceRelation(WVector4_10,WVector4_6,k)-getMultiplicationRecurrenceRelation(WVector4_9,WVector4_5,k);   // Vy_v
//        WVector4_19(k) = getMultiplicationRecurrenceRelation(WVector4_9,WVector4_15,k)-getMultiplicationRecurrenceRelation(XMatrix.row(6),WVector4_7,k)+getMultiplicationRecurrenceRelation(WVector4_10,WVector4_16,k); // Vz_v
//        WVector4_20(k) = getMultiplicationRecurrenceRelation(WVector4_17,WVector4_17,k)+getMultiplicationRecurrenceRelation(WVector4_18,WVector4_18,k);

////        if (WVector4_17 == Eigen::VectorXd::Zero(k+1)){
////            WVector4_21(k) = fabs(WVector4_18(k));
//////            std::cout<<"right1 k = "<<k<<std::endl;
//////            std::cout<<"w4,18 = "<<WVector4_18(k)<<std::endl;
//////            std::cout<<"fabs(0.000347805161450851) = "<<fabs(0.000347805161450851)<<std::endl;
//////            std::cout<<"abs(2) = "<<abs(2)<<std::endl;
//////            std::cout<<"fabs(2.5) = "<<fabs(2.5)<<std::endl;
//////            std::cout<<"fabs(w4,18) = "<<fabs(WVector4_18(k))<<std::endl;
//////            std::cout<<"w4,21 = "<<WVector4_21(k)<<std::endl;
////        }
////        else if (WVector4_18 == Eigen::VectorXd::Zero(k+1)){
////            WVector4_21(k) = fabs(WVector4_17(k));
////            std::cout<<"right2 k = "<<k<<std::endl;
////        }
////        else{
//        WVector4_21(k) = getPowerRecurrenceRelation(WVector4_20,WVector4_21,0.5,k);
////        }

//        WVector4_11(k) = WVector4_20(k)+getMultiplicationRecurrenceRelation(WVector4_19,WVector4_19,k);
//        WVector4_12(k) = getPowerRecurrenceRelation(WVector4_11,WVector4_12,0.5,k);    // Ground velocity

//        WVector4_22(k) = getDivisionRecurrenceRelation(WVector4_18,WVector4_21,WVector4_22,k);  // sin(chi)
//        WVector4_23(k) = getDivisionRecurrenceRelation(WVector4_17,WVector4_21,WVector4_23,k);  // cos(chi)

////        if (k == 1 && -WVector4_19(1)-WVector4_12(1) <= 1e-15){
////            WVector4_24(k) = 0.0;
////        }
////        else{
//        WVector4_24(k) = getDivisionRecurrenceRelation(-WVector4_19,WVector4_12,WVector4_24,k); // sin(gamma)
////        WVector4_24(k) = 0.0;
////        }
//        WVector4_25(k) = getDivisionRecurrenceRelation(WVector4_21,WVector4_12,WVector4_25,k);  // cos(gamma)
        WVector4_26(k) = getCosineRecurrenceRelation(thrustAzimuthVector,WVector4_28,k);
        WVector4_27(k) = getCosineRecurrenceRelation(thrustElevationVector,WVector4_29,k);
        WVector4_28(k) = getSineRecurrenceRelation(thrustAzimuthVector,WVector4_26,k);
        WVector4_29(k) = getSineRecurrenceRelation(thrustElevationVector,WVector4_27,k);
        WVector4_30(k) = getMultiplicationRecurrenceRelation(WVector4_26,WVector4_27,k);

        WVector4_31(k) = getMultiplicationRecurrenceRelation(WVector4_27,WVector4_28,k);
        WVector4_32(k) = getPowerRecurrenceRelation(XMatrix.row(7),WVector4_32,-1.0,k);
        WVector4_33(k) = Thrust_*WVector4_32(k);
//        WVector4_34(k) = getMultiplicationRecurrenceRelation(WVector4_33,WVector4_30,k);
        WVector4_34(k) = WVector4_33(k)*WVector4_27(0)*WVector4_26(0);

        WVector27_1(k) = WVector4_3(k); // Height derivative is just radius
        WVector27_2(k) = getMultiplicationRecurrenceRelation(WVector27_1,WVector27_1,k);
        WVector27_3(k) = getPowerRecurrenceRelation(WVector27_1,WVector27_3,3.0,k);
        WVector27_4(k) = getPowerRecurrenceRelation(WVector27_1,WVector27_4,4.0,k);
        WVector27_5(k) = getPowerRecurrenceRelation(WVector27_1,WVector27_5,5.0,k);
        WVector27_6(k) = getPowerRecurrenceRelation(WVector27_1,WVector27_6,6.0,k);
        WVector27_7(k) = getPowerRecurrenceRelation(WVector27_1,WVector27_7,7.0,k);
        WVector27_8(k) = getPowerRecurrenceRelation(WVector27_1,WVector27_8,8.0,k);
        WVector27_9(k) = getPowerRecurrenceRelation(WVector27_1,WVector27_9,9.0,k);
        WVector27_10(k) = getPowerRecurrenceRelation(WVector27_1,WVector27_10,10.0,k);

        WVector27_11(k) = densityPolyCoefficients(10)*WVector27_10(k)+
                densityPolyCoefficients(9)*WVector27_9(k)+
                densityPolyCoefficients(8)*WVector27_8(k)+
                densityPolyCoefficients(7)*WVector27_7(k)+
                densityPolyCoefficients(6)*WVector27_6(k)+
                densityPolyCoefficients(5)*WVector27_5(k)+
                densityPolyCoefficients(4)*WVector27_4(k)+
                densityPolyCoefficients(3)*WVector27_3(k)+
                densityPolyCoefficients(2)*WVector27_2(k)+
                densityPolyCoefficients(1)*WVector27_1(k);      // Last constant is gone because of derivative
        WVector27_12(k) = getExponentialRecurrenceRelation(WVector27_11,WVector27_12,k);

        // Determine which section of the temperature curve needs to be used and what the corresponding order is
        // Also, because a computer is less than perfect, a small correction is made to the lower bound of the first section to make sure that the initial altitude is still valid

        if ((temperatureAltitudeRanges(0,0)-0.000000000001) <= WVector27_1(0) && WVector27_1(0) < temperatureAltitudeRanges(0,1)){ // Section 1

        WVector27_13(k) = temperaturePolyCoefficients(0,1)*WVector27_1(k);

        }
        else if (temperatureAltitudeRanges(1,0) <= WVector27_1(0) && WVector27_1(0) < temperatureAltitudeRanges(1,1)){  // Section 2

//        WVector27_13(k) = temperaturePolyCoefficients(1,3)*WVector27_3(k)+temperaturePolyCoefficients(1,2)*WVector27_2(k)+temperaturePolyCoefficients(1,1)*WVector27_1(k);
            WVector27_13(k) = temperaturePolyCoefficients(1,2)*WVector27_2(k)+temperaturePolyCoefficients(1,1)*WVector27_1(k);

        }
        else if (temperatureAltitudeRanges(2,0) <= WVector27_1(0) && WVector27_1(0) < temperatureAltitudeRanges(2,1)){ // Section 3

        WVector27_13(k) = temperaturePolyCoefficients(2,6)*WVector27_6(k)+temperaturePolyCoefficients(2,5)*WVector27_5(k)+temperaturePolyCoefficients(2,4)*WVector27_4(k)+
                temperaturePolyCoefficients(2,3)*WVector27_3(k)+temperaturePolyCoefficients(2,2)*WVector27_2(k)+temperaturePolyCoefficients(2,1)*WVector27_1(k);


        }
        else if (temperatureAltitudeRanges(3,0) <= WVector27_1(0) && WVector27_1(0) < temperatureAltitudeRanges(3,1)){ // Section 4

        WVector27_13(k) = temperaturePolyCoefficients(3,8)*WVector27_8(k)+temperaturePolyCoefficients(3,7)*WVector27_7(k)+temperaturePolyCoefficients(3,6)*WVector27_6(k)+
                temperaturePolyCoefficients(3,5)*WVector27_5(k)+temperaturePolyCoefficients(3,4)*WVector27_4(k)+
                  temperaturePolyCoefficients(3,3)*WVector27_3(k)+temperaturePolyCoefficients(3,2)*WVector27_2(k)+temperaturePolyCoefficients(3,1)*WVector27_1(k);

        }
        else if (temperatureAltitudeRanges(4,0) <= WVector27_1(0)){ // Section 5

        WVector27_13(k) = 0;


        }
        else {

        WVector27_13(k) = temperaturePolyCoefficients(0,1)*WVector27_1(k);


        }; // Temperature

        WVector27_14(k) = getPowerRecurrenceRelation((adiabeticIndex*specificGasConstant*WVector27_13),WVector27_14,0.5,k); // Speed of sound
        WVector27_15(k) = getDivisionRecurrenceRelation(WVector4_12,WVector27_14,WVector27_15,k);   // Mach number

        // Determine which section of the drag coefficient curve needs to be used

        for (int i=0; i < 5+1; i++){

            if (dragCoefficientMachRanges(i,0) <= WVector27_15(0) && WVector27_15(0) < dragCoefficientMachRanges(i,1)){

                sectionCD = i;


            }


        };

        WVector27_16(k) = dragCoefficientPolyCoefficients(sectionCD,1)*WVector27_15(k); // C_D value derivative
        WVector27_17(k) = getMultiplicationRecurrenceRelation(WVector4_12,WVector4_12,k);
        WVector27_18(k) = getMultiplicationRecurrenceRelation(WVector27_17,WVector27_16,k);
        WVector27_19(k) = 0.5*referenceArea*getMultiplicationRecurrenceRelation(WVector27_18,WVector27_12,k); // Drag

        WVector4_35(k) = getDivisionRecurrenceRelation(WVector27_19,XMatrix.row(7),WVector4_35,k);
        WVector4_36(k) = WVector4_34(k)-WVector4_35(k);
//        WVector4_36(k) = 0.0; // Test

//        WVector4_37(k) = getMultiplicationRecurrenceRelation(WVector4_33,WVector4_31,k); // Old
        WVector4_37(k) = WVector4_33(k)*WVector4_27(0)*WVector4_28(0); // New
//        WVector4_37(k) = 0.0; // Test
//        WVector4_38(k) = getMultiplicationRecurrenceRelation(WVector4_33,WVector4_29,k); // Old
        WVector4_38(k) = WVector4_33(k)*WVector4_29(0); // New
//        WVector4_38(k) = 0.0; // Test

        WVector4_39(k) = getDivisionRecurrenceRelation((-standardGravitationalParameter*XMatrix.row(1)),XMatrix.row(9),WVector4_39,k);
//        WVector4_40(k) = -getMultiplicationRecurrenceRelation(WVector4_7,WVector4_23,k);

//        WVector4_41(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_24,k);
//        WVector4_42(k) = -getMultiplicationRecurrenceRelation(WVector4_5,WVector4_22,k);
//        WVector4_43(k) = -getMultiplicationRecurrenceRelation(WVector4_5,WVector4_23,k);
//        WVector4_44(k) = -getMultiplicationRecurrenceRelation(WVector4_8,WVector4_25,k);
//        WVector4_45(k) = getMultiplicationRecurrenceRelation(WVector4_40,WVector4_25,k);
//        WVector4_46(k) = getMultiplicationRecurrenceRelation(WVector4_42,WVector4_25,k);
//        WVector4_47(k) = -getMultiplicationRecurrenceRelation(WVector4_13,WVector4_22,k);
//        WVector4_48(k) = getMultiplicationRecurrenceRelation(WVector4_40,WVector4_24,k);
//        WVector4_49(k) = getMultiplicationRecurrenceRelation(WVector4_42,WVector4_24,k);
//        WVector4_50(k) = getMultiplicationRecurrenceRelation(WVector4_6,(WVector4_45+WVector4_41),k)+WVector4_46(k);

//        WVector4_51(k) = getMultiplicationRecurrenceRelation(WVector4_6,(WVector4_48+WVector4_44),k)+WVector4_49(k);
//        WVector4_52(k) = WVector4_39(k)+getMultiplicationRecurrenceRelation(WVector4_36,WVector4_50,k)+getMultiplicationRecurrenceRelation(WVector4_37,(WVector4_47+WVector4_43),k)-
//                getMultiplicationRecurrenceRelation(WVector4_38,WVector4_51,k);

        /// Debug ///
        if (k == 1 || k==2){
//            std::cout<<"k = "<<k<<std::endl;
//        std::cout<<"WVector4_52(k) = "<<WVector4_52(k)<<std::endl;
//        std::cout<<"WVector4_39(k) = "<<WVector4_39(k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_36,WVector4_50,k) = "<<getMultiplicationRecurrenceRelation(WVector4_36,WVector4_50,k)<<std::endl;
////        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_37,(WVector4_47+WVector4_43),k) = "<<getMultiplicationRecurrenceRelation(WVector4_37,(WVector4_47+WVector4_43),k)<<std::endl;
////        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_38,WVector4_51,k) = "<<getMultiplicationRecurrenceRelation(WVector4_38,WVector4_51,k)<<std::endl;
//        std::cout<<"WVector4_36(k) = "<<WVector4_36(k)<<std::endl;
//        std::cout<<"WVector4_50(k) = "<<WVector4_50(k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_6,(WVector4_45+WVector4_41),k) = "<<getMultiplicationRecurrenceRelation(WVector4_6,(WVector4_45+WVector4_41),k)<<std::endl;
//        std::cout<<"WVector4_46(k) = "<<WVector4_46(k)<<std::endl;
//        std::cout<<"WVector4_42(k) = "<<WVector4_42(k)<<std::endl;
//        std::cout<<"WVector4_25(k) = "<<WVector4_25(k)<<std::endl;
//        std::cout<<"WVector4_22(k) = "<<WVector4_22(k)<<std::endl;
//        std::cout<<"WVector4_18(k) = "<<WVector4_18(k)<<std::endl;
//        std::cout<<"WVector4_21(k) = "<<WVector4_21(k)<<std::endl;
//        std::cout<<"WVector4_10(k) = "<<WVector4_10(k)<<std::endl;
//        std::cout<<"WVector4_9(k) = "<<WVector4_9(k)<<std::endl;





//        std::cout<<"WVector4_6(k) = "<<WVector4_6(k)<<std::endl;
//        std::cout<<"WVector4_45(k) = "<<WVector4_45(k)<<std::endl;
//        std::cout<<"WVector4_41(k) = "<<WVector4_41(k)<<std::endl;
//        std::cout<<"WVector4_8(k) = "<<WVector4_8(k)<<std::endl;
//        std::cout<<"WVector4_24(k) = "<<WVector4_24(k)<<std::endl;
//        std::cout<<"WVector4_19(k) = "<<WVector4_19(k)<<std::endl;
//        std::cout<<"WVector4_12(k) = "<<WVector4_12(k)<<std::endl;





        }
//        std::cout<<"WVector4_39(k) = "<<WVector4_39(k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_36,WVector4_50,k) = "<<getMultiplicationRecurrenceRelation(WVector4_36,WVector4_50,k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_37,(WVector4_47+WVector4_43),k) = "<<getMultiplicationRecurrenceRelation(WVector4_37,(WVector4_47+WVector4_43),k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(-WVector4_38,WVector4_51,k) = "<<getMultiplicationRecurrenceRelation(-WVector4_38,WVector4_51,k)<<std::endl;
        /// Debug ///


//        UMatrix(4,k) = WVector4_52(k);
        UMatrix(4,k) = WVector4_39(k)+getMultiplicationRecurrenceRelation(WVector4_36,WVector4_13,k)+getMultiplicationRecurrenceRelation(WVector4_37,WVector4_21,k)-getMultiplicationRecurrenceRelation(WVector4_38,WVector4_24,k); // U4



        // 5

        WVector5_1(k) = getDivisionRecurrenceRelation((-standardGravitationalParameter*XMatrix.row(2)),XMatrix.row(9),WVector5_1,k);
//        WVector5_2(k) = getMultiplicationRecurrenceRelation(WVector4_6,WVector4_22,k);
//        WVector5_3(k) = getMultiplicationRecurrenceRelation(WVector4_5,(WVector4_45+WVector4_41),k)+getMultiplicationRecurrenceRelation(WVector5_2,WVector4_25,k);
//        WVector5_4(k) = -getMultiplicationRecurrenceRelation(WVector4_14,WVector4_22,k)+getMultiplicationRecurrenceRelation(WVector4_6,WVector4_23,k);
//        WVector5_5(k) = getMultiplicationRecurrenceRelation(WVector4_5,(WVector4_48+WVector4_44),k)+getMultiplicationRecurrenceRelation(WVector5_2,WVector4_24,k);
//        WVector5_6(k) = WVector5_1(k)+getMultiplicationRecurrenceRelation(WVector4_36,WVector5_3,k)+getMultiplicationRecurrenceRelation(WVector4_37,WVector5_4,k)-getMultiplicationRecurrenceRelation(WVector4_38,WVector5_5,k);

        /// Debug ///
        if (k==1 || k==2){
//        std::cout<<"w5,6(k) = "<<WVector5_6(k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_36,WVector5_3,k) = "<<getMultiplicationRecurrenceRelation(WVector4_36,WVector5_3,k)<<std::endl;
//        std::cout<<"w4,36(k) = "<<WVector4_36(k)<<std::endl;
//        std::cout<<"w5,3(k) = "<<WVector5_3(k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_5,(WVector4_45+WVector4_41),k) = "<<getMultiplicationRecurrenceRelation(WVector4_5,(WVector4_45+WVector4_41),k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector5_2,WVector4_25,k) = "<<getMultiplicationRecurrenceRelation(WVector5_2,WVector4_25,k)<<std::endl;
//        std::cout<<"x2 = "<<XMatrix.row(2)<<std::endl;
//        std::cout<<"w4,4 = "<<WVector4_4<<std::endl;
//        std::cout<<"w4,4(1) * w4,5(0) =  "<<WVector4_4(1)*WVector4_5(0)<<std::endl;
//        std::cout<<"x2(1) - w4,4(1) * w4,5(0) =  "<<XMatrix(2,1)-WVector4_4(1)*WVector4_5(0)<<std::endl;
//        std::cout<<"w4,5 = "<<WVector4_5<<std::endl;
//        std::cout<<"w4,6 = "<<WVector4_6<<std::endl;
//        std::cout<<"w4_8 = "<<WVector4_8<<std::endl;
//        std::cout<<"w4_9 = "<<WVector4_9<<std::endl;
//        std::cout<<"w4_10 = "<<WVector4_10<<std::endl;
//        std::cout<<"w4_15 = "<<WVector4_15<<std::endl;
//        std::cout<<"w4_16 = "<<WVector4_16<<std::endl;
//        std::cout<<"w5,2 = "<<WVector5_2<<std::endl;
//        std::cout<<"w4_17 = "<<WVector4_17<<std::endl;
//        std::cout<<"w4_18 = "<<WVector4_18<<std::endl;
//        std::cout<<"w4_21 = "<<WVector4_21<<std::endl;
//        std::cout<<"w4_22 = "<<WVector4_22<<std::endl;
//        std::cout<<"w4_24 = "<<WVector4_24<<std::endl;
//        std::cout<<"w4_25 = "<<WVector4_25<<std::endl;
//        std::cout<<"w4_12 = "<<WVector4_12<<std::endl;
//        std::cout<<"-w4_19 = "<<-WVector4_19<<std::endl;

//        std::cout<<"w4_41 = "<<WVector4_41<<std::endl;

}

        /// Debug ///


//        UMatrix(5,k) = WVector5_6(k);

        UMatrix(5,k) = WVector5_1(k)+getMultiplicationRecurrenceRelation(WVector4_36,WVector4_14,k)+getMultiplicationRecurrenceRelation(WVector4_37,WVector4_22,k)-getMultiplicationRecurrenceRelation(WVector4_38,WVector4_25,k); // U5

//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_36,WVector4_14,k) = "<<getMultiplicationRecurrenceRelation(WVector4_36,WVector4_14,k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_37,WVector4_22,k) = "<<getMultiplicationRecurrenceRelation(WVector4_37,WVector4_22,k)<<std::endl;
//        std::cout<<"-getMultiplicationRecurrenceRelation(WVector4_38,WVector4_25,k) = "<<-getMultiplicationRecurrenceRelation(WVector4_38,WVector4_25,k)<<std::endl;



        // 6

        WVector6_1(k) = getDivisionRecurrenceRelation((-standardGravitationalParameter*XMatrix.row(3)),XMatrix.row(9),WVector6_1,k);
//        WVector6_2(k) = getMultiplicationRecurrenceRelation(WVector4_7,WVector4_24,k);
//        WVector6_3(k) = getMultiplicationRecurrenceRelation(WVector4_8,WVector4_22,k);
//        WVector6_4(k) = -getMultiplicationRecurrenceRelation(WVector4_7,WVector4_25,k);
//        WVector6_5(k) = -getMultiplicationRecurrenceRelation(WVector4_44,WVector4_23,k)+WVector6_2(k);
//        WVector6_6(k) = getMultiplicationRecurrenceRelation(WVector4_41,WVector4_23,k)+WVector6_4(k);
//        WVector6_7(k) = WVector6_1(k)+getMultiplicationRecurrenceRelation(WVector4_36,WVector6_5,k)-getMultiplicationRecurrenceRelation(WVector4_37,WVector6_3,k)-getMultiplicationRecurrenceRelation(WVector4_38,WVector6_6,k);


//        UMatrix(6,k) = WVector6_7(k);
        UMatrix(6,k) = WVector6_1(k)+getMultiplicationRecurrenceRelation(WVector4_36,WVector4_15,k)+getMultiplicationRecurrenceRelation(WVector4_37,WVector4_23,k)-getMultiplicationRecurrenceRelation(WVector4_38,WVector4_40,k);

//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_36,WVector4_15,k) = "<<getMultiplicationRecurrenceRelation(WVector4_36,WVector4_15,k)<<std::endl;
//        std::cout<<"getMultiplicationRecurrenceRelation(WVector4_37,WVector4_23,k) = "<<getMultiplicationRecurrenceRelation(WVector4_37,WVector4_23,k)<<std::endl;
//        std::cout<<"-getMultiplicationRecurrenceRelation(WVector4_38,WVector4_40,k) = "<<-getMultiplicationRecurrenceRelation(WVector4_38,WVector4_40,k)<<std::endl;
//        std::cout<<" WVector4_36(k) = "<<WVector4_36(k)<<std::endl;
//        std::cout<<" WVector4_15(k) = "<<WVector4_15(k)<<std::endl;
//        std::cout<<" WVector4_37(k) = "<<WVector4_37(k)<<std::endl;
//        std::cout<<" WVector4_23(k) = "<<WVector4_23(k)<<std::endl;
//        std::cout<<" WVector4_38(k) = "<<WVector4_38(k)<<std::endl;
//        std::cout<<" WVector4_40(k) = "<<WVector4_40(k)<<std::endl;

//        std::cout<<" "<<std::endl;
//        std::cout<<" "<<std::endl;




        // 7

        UMatrix(7,k) = 0;

        // 8

        WVector8_1(k) = getMultiplicationRecurrenceRelation(XMatrix.row(1),XMatrix.row(4),k);
        WVector8_2(k) = getMultiplicationRecurrenceRelation(XMatrix.row(2),XMatrix.row(5),k);
        WVector8_3(k) = getMultiplicationRecurrenceRelation(XMatrix.row(3),XMatrix.row(6),k);

        UMatrix(8,k) = 2.0*(WVector8_1(k)+WVector8_2(k)+WVector8_3(k));

        // 9

        W9IntermediateVector(k) = getMultiplicationRecurrenceRelation(XMatrix.row(9),UMatrix.row(8),k);

        WVector9(k) = getDivisionRecurrenceRelation(W9IntermediateVector,XMatrix.row(8),WVector9,k);

        UMatrix(9,k) = 1.5*WVector9(k);





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

    stateTaylorCoefficients.row(0) = WVector27_15;
    stateTaylorCoefficients(0,maxOrder) = 0.0;  // Last element is not used so set to zero

    stateTaylorCoefficients.row(1) = XMatrix.row(1);
    stateTaylorCoefficients.row(2) = XMatrix.row(2);
    stateTaylorCoefficients.row(3) = XMatrix.row(3);
    stateTaylorCoefficients.row(4) = XMatrix.row(4);
    stateTaylorCoefficients.row(5) = XMatrix.row(5);
    stateTaylorCoefficients.row(6) = XMatrix.row(6);
    stateTaylorCoefficients.row(7) = XMatrix.row(7);


    /// Debug ///

//    std::cout<<"WVector4_1 = "<<WVector4_1<<std::endl;     // W4,1
//    std::cout<<"WVector4_2 = "<<WVector4_2<<std::endl;     // W4,2
////    std::cout<<"XVector.row(3) = "<<XMatrix.row(3)<<std::endl; // x3
//    std::cout<<"WVector4_3 = "<<WVector4_3<<std::endl;     // W4,3
//    std::cout<<"WVector4_4 = "<<WVector4_4<<std::endl;     // W4,4
//    std::cout<<"WVector4_5 = "<<WVector4_5<<std::endl;     // W4,5
//    std::cout<<"WVector4_6 = "<<WVector4_6<<std::endl;     // W4,6
//    std::cout<<"WVector4_7 = "<<WVector4_7<<std::endl;     // W4,7
//    std::cout<<"WVector4_8 = "<<WVector4_8<<std::endl;     // W4,8
//    std::cout<<"WVector4_9 = "<<WVector4_9<<std::endl;     // W4,9
//    std::cout<<"WVector4_10 = "<<WVector4_10<<std::endl;     // W4,10

//    std::cout<<"WVector4_11 = "<<WVector4_11<<std::endl;     // W4,11
//    std::cout<<"WVector4_12 = "<<WVector4_12<<std::endl;     // W4,12
//    std::cout<<"WVector4_13 = "<<WVector4_13<<std::endl;     // W4,13
//    std::cout<<"WVector4_14 = "<<WVector4_14<<std::endl;     // W4,14
//    std::cout<<"WVector4_15 = "<<WVector4_15<<std::endl;     // W4,15
//    std::cout<<"WVector4_16 = "<<WVector4_16<<std::endl;     // W4,16
//    std::cout<<"WVector4_17 = "<<WVector4_17<<std::endl;     // W4,17
//    std::cout<<"WVector4_18 = "<<WVector4_18<<std::endl;     // W4,18
//    std::cout<<"WVector4_19 = "<<WVector4_19<<std::endl;     // W4,19
//    std::cout<<"WVector4_20 = "<<WVector4_20<<std::endl;     // W4,20

//    std::cout<<"WVector4_21 = "<<WVector4_21<<std::endl;     // W4,21
//    std::cout<<"WVector4_22 = "<<WVector4_22<<std::endl;     // W4,22
//    std::cout<<"WVector4_23 = "<<WVector4_23<<std::endl;     // W4,23
//    std::cout<<"WVector4_24 = "<<WVector4_24<<std::endl;     // W4,24
//    std::cout<<"WVector4_25 = "<<WVector4_25<<std::endl;     // W4,25
//    std::cout<<"WVector4_26 = "<<WVector4_26<<std::endl;     // W4,26
//    std::cout<<"WVector4_27 = "<<WVector4_27<<std::endl;     // W4,27
//    std::cout<<"WVector4_28 = "<<WVector4_28<<std::endl;     // W4,28
//    std::cout<<"WVector4_29 = "<<WVector4_29<<std::endl;     // W4,29
//    std::cout<<"WVector4_30 = "<<WVector4_30<<std::endl;     // W4,30

//    std::cout<<"WVector4_31 = "<<WVector4_31<<std::endl;     // W4,31
//    std::cout<<"WVector4_32 = "<<WVector4_32<<std::endl;     // W4,32
//    std::cout<<"WVector4_33 = "<<WVector4_33<<std::endl;     // W4,33
//    std::cout<<"WVector4_34 = "<<WVector4_34<<std::endl;     // W4,34
//    std::cout<<"WVector4_35 = "<<WVector4_35<<std::endl;     // W4,35
//    std::cout<<"WVector4_36 = "<<WVector4_36<<std::endl;     // W4,36
    std::cout<<"WVector4_37 = "<<WVector4_37<<std::endl;     // W4,37
    std::cout<<"WVector4_38 = "<<WVector4_38<<std::endl;     // W4,38
//    std::cout<<"WVector4_39 = "<<WVector4_39<<std::endl;     // W4,39
//    std::cout<<"WVector4_40 = "<<WVector4_40<<std::endl;     // W4,40

//    std::cout<<"WVector4_41 = "<<WVector4_41<<std::endl;     // W4,41
//    std::cout<<"WVector4_42 = "<<WVector4_42<<std::endl;     // W4,42
//    std::cout<<"WVector4_43 = "<<WVector4_43<<std::endl;     // W4,43
//    std::cout<<"WVector4_44 = "<<WVector4_44<<std::endl;     // W4,44
//    std::cout<<"WVector4_45 = "<<WVector4_45<<std::endl;     // W4,45
//    std::cout<<"WVector4_46 = "<<WVector4_46<<std::endl;     // W4,46
//    std::cout<<"WVector4_47 = "<<WVector4_47<<std::endl;     // W4,47
//    std::cout<<"WVector4_48 = "<<WVector4_48<<std::endl;     // W4,48
//    std::cout<<"WVector4_49 = "<<WVector4_49<<std::endl;     // W4,49
//    std::cout<<"WVector4_50 = "<<WVector4_50<<std::endl;     // W4,50

//    std::cout<<"WVector4_51 = "<<WVector4_51<<std::endl;     // W4,51
//    std::cout<<"WVector4_52 = "<<WVector4_52<<std::endl;     // W4,52

//    std::cout<<"WVector5_1 = "<<WVector5_1<<std::endl;     // W5,1
//    std::cout<<"WVector5_2 = "<<WVector5_2<<std::endl;     // W5,2
//    std::cout<<"WVector5_3 = "<<WVector5_3<<std::endl;     // W5,3
//    std::cout<<"WVector5_4 = "<<WVector5_4<<std::endl;     // W5,4
//    std::cout<<"WVector5_5 = "<<WVector5_5<<std::endl;     // W5,5
//    std::cout<<"WVector5_6 = "<<WVector5_6<<std::endl;     // W5,6

//    std::cout<<"WVector6_1 = "<<WVector6_1<<std::endl;     // W6,1
//    std::cout<<"WVector6_2 = "<<WVector6_2<<std::endl;     // W6,2
//    std::cout<<"WVector6_3 = "<<WVector6_3<<std::endl;     // W6,3
//    std::cout<<"WVector6_4 = "<<WVector6_4<<std::endl;     // W6,4
//    std::cout<<"WVector6_5 = "<<WVector6_5<<std::endl;     // W6,5
//    std::cout<<"WVector6_6 = "<<WVector6_6<<std::endl;     // W6,6
//    std::cout<<"WVector6_7 = "<<WVector6_7<<std::endl;     // W6,7



//    std::cout<<"WVector8_1 = "<<WVector8_1<<std::endl;     // W8,1
//    std::cout<<"WVector8_2 = "<<WVector8_2<<std::endl;     // W8,2
//    std::cout<<"WVector8_3 = "<<WVector8_3<<std::endl;     // W8,3


//    std::cout<<"WVector9 = "<<WVector9<<std::endl;     // W9




//    std::cout<<"WVector27_1 = "<<WVector27_1<<std::endl;     // W27,1
//    std::cout<<"WVector27_2 = "<<WVector27_2<<std::endl;     // W27,2
//    std::cout<<"WVector27_3 = "<<WVector27_3<<std::endl;     // W27,3
//    std::cout<<"WVector27_4 = "<<WVector27_4<<std::endl;     // W27,4
//    std::cout<<"WVector27_5 = "<<WVector27_5<<std::endl;     // W27,5
//    std::cout<<"WVector27_6 = "<<WVector27_6<<std::endl;     // W27,6
//    std::cout<<"WVector27_7 = "<<WVector27_7<<std::endl;     // W27,7
//    std::cout<<"WVector27_8 = "<<WVector27_8<<std::endl;     // W27,8
//    std::cout<<"WVector27_9 = "<<WVector27_9<<std::endl;     // W27,9
//    std::cout<<"WVector27_10 = "<<WVector27_10<<std::endl;     // W27,10

//    std::cout<<"WVector27_11 = "<<WVector27_11<<std::endl;     // W27,11
//    std::cout<<"WVector27_12 = "<<WVector27_12<<std::endl;     // W27,12
//    std::cout<<"WVector27_13 = "<<WVector27_13<<std::endl;     // W27,13
//    std::cout<<"WVector27_14 = "<<WVector27_14<<std::endl;     // W27,14
//    std::cout<<"WVector27_15 = "<<WVector27_15<<std::endl;     // W27,15
//    std::cout<<"WVector27_16 = "<<WVector27_16<<std::endl;     // W27,16
//    std::cout<<"WVector27_17 = "<<WVector27_17<<std::endl;     // W27,17
//    std::cout<<"WVector27_18 = "<<WVector27_18<<std::endl;     // W27,18
//    std::cout<<"WVector27_19 = "<<WVector27_19<<std::endl;     // W27,19

//Eigen::VectorXd diff = Eigen::VectorXd::Zero(maxOrder);
//Eigen::VectorXd fraction = Eigen::VectorXd::Zero(maxOrder);
//for (int i = 0; i < maxOrder; i++){
//       diff(i) = XMatrix(3,i)-WVector4_3(i);
//       fraction(i) = XMatrix(3,i)/WVector4_3(i)-1;
//}
//std::cout<<"x3-w4,3 = "<<diff<<std::endl;
//std::cout<<"x3/w4,3-1 = "<<fraction<<std::endl;

//    std::cout<<"XMatrix.row(6) = "<<XMatrix.row(6)<<std::endl;


//    std::cout<<"XMatrix = "<<XMatrix<<std::endl;
//    std::cout<<"UMatrix = "<<UMatrix<<std::endl;






    return stateTaylorCoefficients;

} // end of function




