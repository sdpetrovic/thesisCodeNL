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
 *      160504    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */




#ifndef STEPSIZE_H
#define STEPSIZE_H


#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>

#include <tudat/Tudat/Mathematics/BasicMathematics/linearAlgebra.h>
#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>


class StepSize

        /* This class will provide the current step-size, a manner to manually update/set the step-size, local error tolerance and step multiplication factor.
         * It will also provide the truncation error estimate and the truncation error per state variable.
         *
         * The equations provided in this class are all taken from the following paper:
         *
         * "High Speed Solution of Spacecraft Trajectory Problems Using Taylor Series Integration"
         * by James R. Scott and Michael C. Martini
         * AIAA/AAS Astrodynamics Specialist Conference and Exhibit
         * 18-21 August 2008, Honolulu, Hawaii
         * AIAA 2008-6957
         * ///////////////////////////////////////////////////////////////////////////////////////
         *
         * In the original formulation the different variables are represented as:
         *
         * - h [s]                      step-size
         * - h_next [s]                 next step-size
         * - eta                        step multiplication factor
         * - tau                        local error tolerance
         * - M                          order of the maximum truncation error estimate
         * - e_max                      maximum truncation error estimate
         * - T_n,K                      state variable truncation error estimate
         *
         */


{
public:


    // Initializing the StepSize class with de default setting for the step-size set to 1 sec

    StepSize(const double currentStepSize_ = 1){

        currentStepSize = currentStepSize_;             // The current step-size
        stepMultiplicationFactor = 0.9;                 // Setting the step multiplication factor to the default value
        localErrorTolerance = 1E-2;                     // Setting the local error tolerance to the default value


    } // End of constructor


    /// All the set functions for the three variables mentioned in the class construction ///

    // Current step-size
    void setCurrentStepSize( const double updatedStepSize  )
    {
        currentStepSize = updatedStepSize;
    }

    // Step multiplication factor
    void setStepMultiplicationFactor( const double updatedStepMultiplicationFactor ){

        stepMultiplicationFactor = updatedStepMultiplicationFactor;
    }

    // Local error tolerance
    void setLocalErrorTolerance( const double updatedLocalErrorTolerance ){

        localErrorTolerance = updatedLocalErrorTolerance;
    }


    /// Determine the estimate of the truncation error ///

    void determineTruncationErrorEstimates(const tudat::basic_mathematics::Vector7d& penultimateCoefficients,
                                                                               const tudat::basic_mathematics::Vector7d& lastCoefficients, const int maxOrder){

        truncationErrorEstimates(0) = abs(penultimateCoefficients(0))*pow(currentStepSize,(maxOrder-1))+abs(lastCoefficients(0))*pow(currentStepSize,maxOrder);     // Position in x-direction
        truncationErrorEstimates(1) = abs(penultimateCoefficients(1))*pow(currentStepSize,(maxOrder-1))+abs(lastCoefficients(1))*pow(currentStepSize,maxOrder);     // Position in y-direction
        truncationErrorEstimates(2) = abs(penultimateCoefficients(2))*pow(currentStepSize,(maxOrder-1))+abs(lastCoefficients(2))*pow(currentStepSize,maxOrder);     // Position in z-direction
        truncationErrorEstimates(3) = abs(penultimateCoefficients(3))*pow(currentStepSize,(maxOrder-1))+abs(lastCoefficients(3))*pow(currentStepSize,maxOrder);     // Velocity in x-direction
        truncationErrorEstimates(4) = abs(penultimateCoefficients(4))*pow(currentStepSize,(maxOrder-1))+abs(lastCoefficients(4))*pow(currentStepSize,maxOrder);     // Velocity in y-direction
        truncationErrorEstimates(5) = abs(penultimateCoefficients(5))*pow(currentStepSize,(maxOrder-1))+abs(lastCoefficients(5))*pow(currentStepSize,maxOrder);     // Velocity in z-direction
        truncationErrorEstimates(6) = abs(penultimateCoefficients(6))*pow(currentStepSize,(maxOrder-1))+abs(lastCoefficients(6))*pow(currentStepSize,maxOrder);     // Mass

    }


    /// Determine the estimate of the maximum truncation error ///




    void determineMaximumTruncationErrorEstimate(const tudat::basic_mathematics::Vector7d& truncationErrorEstimates_){

        double max; // Define the maximum holding variable

        for (int i = 0; i<truncationErrorEstimates_.size()-1;i++){

            if (truncationErrorEstimates_(i)>=truncationErrorEstimates_(i+1)){

                max = truncationErrorEstimates_(i);

            } else if (truncationErrorEstimates_(i)<truncationErrorEstimates_(i+1)){

                max = truncationErrorEstimates_(i+1);
            }
        }

        maximumTruncationErrorEstimate = max;

//        maximumTruncationErrorEstimate = max(truncationErrorEstimates_);


    }


    /// Determine the next step-size using the previous step-size///

    void determineNextStepSizeUsingPreviousStepSize(const double maxOrder_){


        // Determine the order of magnitude of the maximum truncation error
         orderMaxTruncErrorEstimate = maxOrder_;


        nextStepSize = stepMultiplicationFactor*currentStepSize*pow((localErrorTolerance/maximumTruncationErrorEstimate),(1/orderMaxTruncErrorEstimate));


/*        /// Debug ///

        std::cout<<"localErrorTolerance/maximumTruncationErrorEstimate = "<<localErrorTolerance/maximumTruncationErrorEstimate<<std::endl;
        std::cout<<"1/orderMaxTruncErrorEstimate = "<<(1/orderMaxTruncErrorEstimate)<<std::endl;
        std::cout<<"pow((localErrorTolerance/maximumTruncationErrorEstimate),(1/orderMaxTruncErrorEstimate)) = "<<pow((localErrorTolerance/maximumTruncationErrorEstimate),(1/orderMaxTruncErrorEstimate))<<std::endl;
//*/
    }


    /// Determine the next step-size using iteration and the local error tolerance directly ///

    void determineNextStepSizeUsingIteration(const tudat::basic_mathematics::Vector7d& penultimateCoefficients_,
                                             const tudat::basic_mathematics::Vector7d& lastCoefficients_, const int maxOrder_){

        double maxOrder = maxOrder_; // Making it a double such that it can be used for computations
        tudat::basic_mathematics::Vector7d allStepSizes;        // Vector containing the step-sizes for each variable

        for (int i = 0; i<penultimateCoefficients_.size();i++){             // Determining the required step-size for every variable

            double prevStepSize = 1;        // Default
            double newStepSize = 0;         // Default

            double differenceIntStepSize = prevStepSize-newStepSize;        // Difference to determine if the while loop needs to keep going



            while (abs(differenceIntStepSize) > 1e-15){                      // Accepting the step-size if the difference is 1e-6

                    newStepSize = exp((1/(maxOrder-1))*log(localErrorTolerance/(abs(penultimateCoefficients_(i))+prevStepSize*abs(lastCoefficients_(i)))));  // Eq. 47 of the reference paper (High Speed Solution of Spacecraft Trajectory Problems using Taylor Series Integration

                    differenceIntStepSize = prevStepSize-newStepSize;           // Determining the difference

//                    std::cout<<"differenceIntstepSize = "<<differenceIntStepSize<<std::endl;

                    prevStepSize = newStepSize;

//                    std::cout<<"newStepSize = "<<newStepSize<<std::endl;
            }

            allStepSizes(i) = prevStepSize;     // Storing the step-size

        }

        double min; // Define the minimum holding variable

        for (int i = 0; i<allStepSizes.size()-1;i++){

            if (allStepSizes(i)>=allStepSizes(i+1)){

                min = allStepSizes(i+1);

            } else if (allStepSizes(i)<allStepSizes(i+1)){

                min = allStepSizes(i);
            }
        }

        std::cout<<"allStepSizes = "<<allStepSizes<<std::endl;
        std::cout<<"min = "<<min<<std::endl;

        nextStepSize = stepMultiplicationFactor*min;

    }


    /// One overall function to update the next step-size (using the above mentioned functions) ///


    // Using the previous step-size
    void updateStepSizeUsingPreviousStepSize(const tudat::basic_mathematics::Vector7d& penultimateCoefficients_,
                        const tudat::basic_mathematics::Vector7d& lastCoefficients_, const int maxOrder_){

        determineTruncationErrorEstimates(penultimateCoefficients_, lastCoefficients_, maxOrder_);
        determineMaximumTruncationErrorEstimate(truncationErrorEstimates);
        determineNextStepSizeUsingPreviousStepSize(maxOrder_);

        // Updating the step-size
        currentStepSize = nextStepSize;
    }


    // Using the iteration method
    void updateStepSizeUsingIteration(const tudat::basic_mathematics::Vector7d& penultimateCoefficients_,
                                     const tudat::basic_mathematics::Vector7d& lastCoefficients_, const int maxOrder_){

        determineNextStepSizeUsingIteration(penultimateCoefficients_, lastCoefficients_, maxOrder_);

        // Updating the step-size
        currentStepSize = nextStepSize;
    }



    /// All return functions ///

    const double getCurrentStepSize() { return currentStepSize;}                     // The current step-size (or h)
    const double getStepMultiplicationFactor() { return stepMultiplicationFactor;}             // The step multiplication factor (or eta)
    const double getLocalErrorTolerance() {return localErrorTolerance;}                 // The local error tolerance (or tau)


    const double getNextStepSize() {return nextStepSize;}                               // The next step-size (or h_next)
    const double getOrderMaxTruncErrorEstimate() {return orderMaxTruncErrorEstimate;}          // The order of the maximum truncation error estimate (or M)
    const double getMaximumTruncationErrorEstimate(){return maximumTruncationErrorEstimate;}      // The maximum truncation error estimate (or e_max)


    const tudat::basic_mathematics::Vector7d getTruncationErrorEstimates(){return truncationErrorEstimates;}         // The state variable truncation error estimates or (T_n,K)




private:

        // Defining the different parameters

    double currentStepSize;                     // The current step-size (or h)
    double stepMultiplicationFactor;             // The step multiplication factor (or eta)
    double localErrorTolerance;                 // The local error tolerance (or tau)


    double nextStepSize;                        // The next step-size (or h_next)
    double orderMaxTruncErrorEstimate;          // The order of the maximum truncation error estimate (or M)
    double maximumTruncationErrorEstimate;      // The maximum truncation error estimate (or e_max)


    tudat::basic_mathematics::Vector7d truncationErrorEstimates;         // The state variable truncation error estimates or (T_n,K)

};

#endif // STEPSIZE_H
