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
 *      121030    K. Kumar          File created.
 *      130225    K. Kumar          Updated include-guard and namespace names; updated Vector6d
 *                                  references to use Tudat definition.
 *      160412    S.D. Petrovic     Modified for thesis application (change name of body class and added variable)
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef STATE_AND_TIME_H
#define STATE_AND_TIME_H

#include <map>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
//#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <tudatApplications/thesisProject/linearAlgebraTypesUpdated.h>



//! Test body class.
/*!
 * This class serves as an example of how a container can be constructed that stores state and
 * time information, which can be used in conjunction with acceleration models, the Cartesian state
 * derivative model, and the composite state derivative model. It should be noted that this class
 * should NOT be used "as is", without consideration for the application at hand. Classes such as
 * this are application-specific, hence unavailable through the Tudat libraries.
 */
class StateAndTime
{
public:

    //! Constructor taking a state (,mass) and a time.
    /*!
     * Constructor taking an input state (mass separate) and time. The input state is used internally to
     * set the current position (taken as a segment of the input state given by the indices
     * (0, 3)), the current velocity (taken as a segment of the input state given by the indices
     * (3, 3)) and the current mass (taken as a segment of the input state given by the indices
     * (6, 1)).
     * \param aState An input state vector.
     * \param aTime An input time (default = 0.0) [s].
     */
    StateAndTime( const tudat::basic_mathematics::Vector7d& aState, const double aTime = 0.0 )
        : currentState( aState ),
          currentPosition( aState.segment( 0, 3 ) ),
          currentVelocity( aState.segment( 3, 3 ) ),
          currentTime( aTime ),
          currentMass( aState(6) )
    { }

    //! Set current time and state.
    /*!
     * Sets the current time, position and current velocity internally based on the input
     * arguments. The current position is taken as a segment of the input state given by the
     * indices (0, 3), the current velocity is taken as a segment of the input state given by
     * the indices (3, 3) and the current mass is taken as a segment of the input state given
     * by the indices (6, 1).
     * \param aTime An input time [s].
     * \param aState An input state vector.
     */
    void setCurrentStateAndTime( const tudat::basic_mathematics::Vector7d& aState, const double aTime  )
    {
        currentTime = aTime;
        currentState = aState;
        currentPosition = aState.segment( 0, 3 );
        currentVelocity = aState.segment( 3, 3 );
        currentMass = aState(6);
    }

    //! Get current state.
    /*!
     * Returns the internally stored current state vector.
     * \return Current state.
     */
    tudat::basic_mathematics::Vector7d getCurrentState( ) { return currentState; }

    //! Get current position.
    /*!
     * Returns the internally stored current position vector.
     * \return Current position.
     */
    Eigen::Vector3d getCurrentPosition( ) { return currentPosition; }

    //! Get current velocity.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity.
     */
    Eigen::Vector3d getCurrentVelocity( ) { return currentVelocity; }

    //! Get current time.
    /*!
     * Returns the internally stored current time.
     * \return Current time.
     */
    double getCurrentMass( ) { return currentMass; }




    //! Get current time.
    /*!
     * Returns the internally stored current time.
     * \return Current time.
     */
    double getCurrentTime( ) { return currentTime; }

protected:

private:

    //! Current state.
    tudat::basic_mathematics::Vector7d currentState;

    //! Current position.
    Eigen::Vector3d currentPosition;

    //! Current position.
    Eigen::Vector3d currentVelocity;

    //! Current time.
    double currentTime;

    //! Current mass.
    double currentMass;


};

//! Typedef for shared-pointer to stateAndTime.
typedef boost::shared_ptr< StateAndTime > StateAndTimePointer;


#endif // STATE_AND_TIME_H
