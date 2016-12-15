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
 *      160427    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef BASICRECURRENCERELATIONS_H
#define BASICRECURRENCERELATIONS_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>

//class basicRecurrenceRelations
//{
//public:
//    basicRecurrenceRelations();
//};



////// Declare functions ///////

/// Multiplication ///
/// \brief getMultiplicationRecurrenceRelation
/// \param F
/// \param G
/// \return
///
///
double getMultiplicationRecurrenceRelation(
       const Eigen::VectorXd& F,
       const Eigen::VectorXd& G,
        const int order);


/// Division ///
/// \brief getDivisionRecurrenceRelation
/// \param F
/// \param G
/// \param Wdiv
/// \return
///
double getDivisionRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& G,
        const Eigen::VectorXd& Wdiv,
        const int order);


/// Power ///
/// \brief getPowerRecurrenceRelation
/// \param F
/// \param Wpow
/// \param power
/// \return
///
double getPowerRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& Wpow,
        const double power,
        const int order);


/// Exponential ///
/// \brief getExponentialRecurrenceRelation
/// \param F
/// \param Wexp
/// \return
///
double getExponentialRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& Wexp,
        const int order);


/// Cosine ///
/// \brief getCosineRecurrenceRelation
/// \param F
/// \param Wsin
/// \return
///
double getCosineRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& Wsin,
        const int order);


/// Sine ///
/// \brief getSineRecurrenceRelation
/// \param F
/// \param Wcos
/// \return
///
double getSineRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& Wcos,
        const int order);





#endif // BASICRECURRENCERELATIONS_H
