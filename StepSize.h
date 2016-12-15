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
 *      160505    S.D. Petrovic     Redefined M and added second step-size determination method. Also moved some of the calculations into private
 *      160531    S.D. Petrovic     Made changes to the iteration method using Bergsma's thesis
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
#include <complex>
#include <Eigen/Core>

#include <tudat/Tudat/Mathematics/BasicMathematics/linearAlgebra.h>
#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>


class StepSize

        /* This class will provide the current step-size, a manner to manually update/set the step-size, local error tolerance and step multiplication factor.
         * It will provide two different methods to determine the next step-size. But only one of these methods should be used at once!
         *
         * ///////////////////////////////////////////////////////////////////////////////////////
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

    StepSize(const double currentStepSize_ = 0.2){

        currentStepSize = currentStepSize_;             // The current step-size
        stepMultiplicationFactor = 0.9;                 // Setting the step multiplication factor to the default value
        localErrorTolerance = 1E-14;                     // Setting the local error tolerance to the default value


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

/*///// Calculations ////////////

    /// Determine the estimate of the truncation error ///

    void determineTruncationErrorEstimates(const tudat::basic_mathematics::Vector7d& penultimateCoefficients,
                                                                               const tudat::basic_mathematics::Vector7d& lastCoefficients, const int maxOrder){         // Eq. 45 of reference paper

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

        max = truncationErrorEstimates_(0);                         // Define the first one to be max initially

        for (int i = 1; i<truncationErrorEstimates_.size();i++){                          //Check for all elements

         if (max<truncationErrorEstimates_(i)){                     // See if current element is larger than the current max

                max = truncationErrorEstimates_(i);                 // If so, set the new max to be the current element value
            }
        }

        maximumTruncationErrorEstimate = max;


    }


    /// Determine the next step-size using the previous step-size///

    void determineNextStepSizeUsingPreviousStepSize(const double maxOrder_){


        // Determine the order of magnitude of the maximum truncation error
         orderMaxTruncErrorEstimate = maxOrder_;        // Apparently this should be the order of the last element of the taylor series, so if I am doing a 20 order taylor series, this should be 20 (?)

        nextStepSize = stepMultiplicationFactor*currentStepSize*pow((localErrorTolerance/maximumTruncationErrorEstimate),(1/orderMaxTruncErrorEstimate)); // Eq. 44 of reference paper


      /// Debug ///

//        std::cout<<"localErrorTolerance/maximumTruncationErrorEstimate = "<<localErrorTolerance/maximumTruncationErrorEstimate<<std::endl;
//        std::cout<<"1/orderMaxTruncErrorEstimate = "<<(1/orderMaxTruncErrorEstimate)<<std::endl;
//        std::cout<<"pow((localErrorTolerance/maximumTruncationErrorEstimate),(1/orderMaxTruncErrorEstimate)) = "<<pow((localErrorTolerance/maximumTruncationErrorEstimate),(1/orderMaxTruncErrorEstimate))<<std::endl;

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

                    prevStepSize = newStepSize;                     // Updating step-size

//                    std::cout<<"newStepSize = "<<newStepSize<<std::endl;
            }

            allStepSizes(i) = prevStepSize;     // Storing the step-size

        }

        double min; // Define the minimum holding variable

        min = allStepSizes(0);                          // Set the minimum value to the first element as default

        for (int i = 1; i<allStepSizes.size();i++){     // Test all other elements

           if (min>allStepSizes(i)){                    // See if the current element is smaller than the current minimum

                min = allStepSizes(i);                  // If so, then the new minimum is set to be the current element value
            }
        }

//        std::cout<<"allStepSizes = "<<allStepSizes<<std::endl;
//        std::cout<<"min = "<<min<<std::endl;

        nextStepSize = stepMultiplicationFactor*min;   // The next step-size then becomes the minimum step-size times the stepMultiplicationFactor

    }
//*/

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

    ///// Calculations ////////////

        /// Determine the estimate of the truncation error ///

        void determineTruncationErrorEstimates(const tudat::basic_mathematics::Vector7d& penultimateCoefficients,
                                                                                   const tudat::basic_mathematics::Vector7d& lastCoefficients, const int maxOrder){         // Eq. 45 of reference paper

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

            max = truncationErrorEstimates_(0);                         // Define the first one to be max initially

            for (int i = 1; i<truncationErrorEstimates_.size();i++){                          //Check for all elements

             if (max<truncationErrorEstimates_(i)){                     // See if current element is larger than the current max

                    max = truncationErrorEstimates_(i);                 // If so, set the new max to be the current element value
                }
            }

            maximumTruncationErrorEstimate = max;


        }


        /// Determine the next step-size using the previous step-size///

        void determineNextStepSizeUsingPreviousStepSize(const double maxOrder_){


            // Determine the order of magnitude of the maximum truncation error
             orderMaxTruncErrorEstimate = maxOrder_;        // Apparently this should be the order of the last element of the taylor series, so if I am doing a 20 order taylor series, this should be 20 (?)

            nextStepSize = stepMultiplicationFactor*currentStepSize*pow((localErrorTolerance/maximumTruncationErrorEstimate),(1/orderMaxTruncErrorEstimate)); // Eq. 44 of reference paper


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
            allStepSizes << 1, 1, 1, 1, 1, 1, 1; // Fill with default value of 1 sec

            for (int i = 0; i<penultimateCoefficients_.size();i++){             // Determining the required step-size for every variable

//                std::cout<<"This is x("<<i<<")"<<std::endl;

                double prevStepSize;
                if (penultimateCoefficients_(i) == 0){
                    prevStepSize = 1;
                }
                else{// From Bergsma's thesis! // This method saves about 1 iteration step...
                    std::complex<double> mycomplex(std::pow(std::complex<double> (localErrorTolerance/penultimateCoefficients_(i)),std::complex<double> (1/(maxOrder-1))));
//                prevStepSize = (pow((localErrorTolerance/penultimateCoefficients_(i)),(1/(maxOrder-1)))); // Didn't work because it resulted in imaginary values sometimes
                    prevStepSize = mycomplex.real();
//                std::cout<<"prevStepSize (using Bergsma's thingy) = "<<prevStepSize<<std::endl;

                // IMPORTANT!!: Sometimes the initial guess from Bergsma gives a number that is maller than the local error tolerance, resulting in an initial stepSize guess
                // that is too small. The actual loop will then provide a larger stepSize guess more proportional to the others. Therefore, if this happens,the initial guess
                // in the form of the prevStepSize is set to 1 (using my initial method) such that the loop will actually take place.

                if (prevStepSize<localErrorTolerance){ // Reset of the prevStepSize in case the initial guess is lower than the localErrorTolerance.
                    prevStepSize = 1.0;
                }

                /// Debug ///
//                std::cout<<"penultimateCoefficients_(i) = "<<penultimateCoefficients_(i)<<std::endl;
//                std::cout<<"(localErrorTolerance/penultimateCoefficients_(i)) = "<<(localErrorTolerance/penultimateCoefficients_(i))<<std::endl;
//                std::cout<<"(1/(maxOrder-1)) = "<<(1/(maxOrder-1))<<std::endl;

                /// Debug ///
                }
                double newStepSize = 0;         // Default

                double differenceIntStepSize = prevStepSize-newStepSize;        // Difference to determine if the while loop needs to keep going
//                std::cout<<"abs(differenceIntStepSize) = "<<abs(differenceIntStepSize)<<std::endl;
                double prevDifferenceIntStepSize = 0;                           // Difference to compare to the current difference to avoid oscillation
                double nNotChanged = 0;             // Set the number of not changed double loops to 0

                if (abs(penultimateCoefficients_(i)) == 0 && abs(lastCoefficients_(i)) == 0){
                    allStepSizes(i) = -log(0);
                }
                else{
                while ((abs(differenceIntStepSize) > 1e-15) && (nNotChanged < 4)){                      // Accepting the step-size if the difference is 1e-15
//                    std::cout<<"nNotchanged = "<<nNotChanged<<std::endl;
//                    std::cout<<"i = "<<i<<std::endl;
//                    std::cout<<"prevStepSize = "<<prevStepSize<<std::endl;
//                    std::cout<<"newStepSize before = "<<newStepSize<<std::endl;

                    if (maxOrder == 1){
                        newStepSize = exp(1.0*log(localErrorTolerance/(abs(penultimateCoefficients_(i))+maxOrder*prevStepSize*abs(lastCoefficients_(i))))); // Changed including maxOrder as by Bergsma (same in the one below)
                    }
                    else{
                        newStepSize = exp((1.0/(maxOrder-1.0))*log(localErrorTolerance/(abs(penultimateCoefficients_(i))+maxOrder*prevStepSize*abs(lastCoefficients_(i)))));  // Eq. 47 of the reference paper (High Speed Solution of Spacecraft Trajectory Problems using Taylor Series Integration
};
                        /// Debug start ///

//                        std::cout<<"newStepSize after = "<<newStepSize<<std::endl;
//                        std::cout<<"1/(maxOrder-1) = "<<1/(maxOrder-1)<<std::endl;
//                        std::cout<<"maxOrder = "<<maxOrder<<std::endl;
//                        std::cout<<"abs(penultimateCoefficients_(i)) = "<<abs(penultimateCoefficients_(i))<<std::endl;
//                        std::cout<<"abs(lastCoefficients_(i)) = "<<abs(lastCoefficients_(i))<<std::endl;
//                        std::cout<<"(abs(penultimateCoefficients_(i))+prevStepSize*abs(lastCoefficients_(i))) = "<<(abs(penultimateCoefficients_(i))+prevStepSize*abs(lastCoefficients_(i)))<<std::endl;
//                        std::cout<<"(localErrorTolerance/(abs(penultimateCoefficients_(i))+prevStepSize*abs(lastCoefficients_(i)))) = "<<(localErrorTolerance/(abs(penultimateCoefficients_(i))+prevStepSize*abs(lastCoefficients_(i))))<<std::endl;
//                        std::cout<<"log(localErrorTolerance/(abs(penultimateCoefficients_(i))+prevStepSize*abs(lastCoefficients_(i)))) = "<<log(localErrorTolerance/(abs(penultimateCoefficients_(i))+prevStepSize*abs(lastCoefficients_(i))))<<std::endl;


                        /// Debug end ///


                        differenceIntStepSize = prevStepSize-newStepSize;           // Determining the difference

                        if (prevDifferenceIntStepSize+differenceIntStepSize == 0){

                            nNotChanged++;

                        }

//                        std::cout<<"prevDifferenceIntStepSize+differenceIntStepSize = "<<prevDifferenceIntStepSize+differenceIntStepSize<<std::endl;
                        prevDifferenceIntStepSize = differenceIntStepSize;          // Update the previous difference

//                        std::cout<<"differenceIntStepSize = "<<differenceIntStepSize<<std::endl;
//                        std::cout<<"prevDifferenceIntStepSize = "<<differenceIntStepSize<<std::endl;


                        prevStepSize = newStepSize;                     // Updating step-size

//                        std::cout<<"newStepSize = "<<newStepSize<<std::endl;
                }   // end while loop

                allStepSizes(i) = prevStepSize;     // Storing the step-size
//                std::cout<<"allStepSizes(i) = "<<allStepSizes(i)<<"     ----------------------------------------------"<<std::endl;

            }       // end else statement
            }   // end for loop

            double min; // Define the minimum holding variable

            min = allStepSizes(0);                          // Set the minimum value to the first element as default

            for (int i = 1; i<allStepSizes.size();i++){     // Test all other elements

               if (min>allStepSizes(i)){                    // See if the current element is smaller than the current minimum

                    min = allStepSizes(i);                  // If so, then the new minimum is set to be the current element value
                }
            }

//            std::cout<<"allStepSizes = "<<allStepSizes<<std::endl;
//            std::cout<<"min = "<<min<<std::endl;

            nextStepSize = stepMultiplicationFactor*min;   // The next step-size then becomes the minimum step-size times the stepMultiplicationFactor
//            std::cout<<"nextStepSize = "<<nextStepSize<<std::endl;

        }

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
