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

#include "basicRecurrenceRelations.h"

//basicRecurrenceRelations::basicRecurrenceRelations()
//{

//}

/// Multiplication ///
double getMultiplicationRecurrenceRelation(
       const Eigen::VectorXd& F,
       const Eigen::VectorXd& G,
        const int order){

    double Wmult_ = 0;                   // Setting the output to zero initially (just in case)

//const int order = F.size()-1;           // Determining the order of the computation, which is the length of the F vector minus one (because the first entry corresponds to the 0th order, etc.)


// Computing the recurrence value
for (int j=0; j < order+1; j++){                    // It goes till the order (till k), and stops as soon as j becomes k+1


     Wmult_ += F(j)*G(order-j);          // Wmult += ... means Wmult = Wmult + ...      // Corrected a problem where I had: F(order)*G(order-j) instead of F(j)*G(order-j)

};

 return Wmult_;


}


/// Division ///
double getDivisionRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& G,
        const Eigen::VectorXd& Wdiv,
        const int order){               // Please note that the previous outcomes should be known

    double Wdiv_ = 0;               // Setting the outcome to zero to avoid previous answers that might be stored somewhere for some reason
    double interSum = 0;           // Defining an intermediate outcome for the summation

    // Avoid singularities

    if (G(0) != 0){

        //    const int order = F.size()-1;               // Determine the current order

            for (int j=1; j<order+1; j++){

            interSum +=G(j)*Wdiv(order-j);              // Perform the summation

            };

            Wdiv_ = (1/G(0))*(F(order)-interSum);       // Compute the current value

    }




    return Wdiv_;

}


/// Power ///
double getPowerRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& Wpow,
        const double power,
        const int order){


    double Wpow_ = 0;               // Setting the outcome to zero to avoid previous answers that might be stored somewhere for some reason
    double interSum = 0;            // Defining an intermediate outcome for the summation


    // Avoid singularities (only compute if the denominator

    if (F(0) != 0){

        //    const int order = F.size()-1;                   // Determine the current order

            for (int j=0; j<order; j++){

                interSum += (order*power-j*(power+1))*F(order-j)*Wpow(j);    // Perform the summation
            };

            Wpow_ = interSum/(order*F(0));          // Compute the current value

        // Same as Bergsma:

//                    for (int j=1; j<order+1; j++){

//                        interSum += ((j/order)*(power+1)-1)*F(j)*Wpow(order-j);    // Perform the summation
//                    };

//                    Wpow_ = interSum/(F(0));          // Compute the current value

    }



    return Wpow_;
}


/// Exponential ///
double getExponentialRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& Wexp,
        const int order){

    double Wexp_ = 0;               // Setting the outcome to zero to avoid previous answers that might be stored somewhere for some reason
    double interSum = 0;            // Defining an intermediate outcome for the summation

//    const int order = F.size()-1;                   // Determine the current order

    for (int j=0; j<order; j++){

        interSum += (order-j)*Wexp(j)*F(order-j);           // Perform the summation
    };

    Wexp_ = interSum/order;                 // Compute the current value

    return Wexp_;

}





/// Cosine ///
double getCosineRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& Wsin,
        const int order){


    double Wcos_ = 0;               // Setting the outcome to zero to avoid previous answers that might be stored somewhere for some reason
    double interSum = 0;            // Defining an intermediate outcome for the summation

//    const int order = F.size()-1;                   // Determine the current order

    for (int j=1; j<order+1; j++){

        interSum +=j*Wsin(order-j)*F(j);        // Perform the summation

    };

    Wcos_ = -interSum/order;                // Compute the current value

    return Wcos_;

}


/// Sine ///
double getSineRecurrenceRelation(
        const Eigen::VectorXd& F,
        const Eigen::VectorXd& Wcos,
        const int order){

    double Wsin_ = 0;               // Setting the outcome to zero to avoid previous answers that might be stored somewhere for some reason
    double interSum = 0;            // Defining an intermediate outcome for the summation

//    const int order = F.size()-1;                   // Determine the current order

    for (int j=1; j<order+1; j++){

        interSum +=j*Wcos(order-j)*F(j);        // Perform the summation

    };

    Wsin_ = interSum/order;                // Compute the current value

    return Wsin_;

}
