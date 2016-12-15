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
 *      160516    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

#include "otherRequiredFunctions.h"


///// Some required functions /////

const long double pi = 4.0*atan(1);


/// deg2rad and rad2deg ///

const double deg2rad(const double deg){

//    const double rad = deg*tudat::mathematical_constants::LONG_PI/180;
    const double rad = deg*pi/180.0;
//    std::cout<<"pi = "<<pi<<std::endl;
//    std::cout<<"PI = "<<tudat::mathematical_constants::LONG_PI<<std::endl;

    return rad;
}

const double rad2deg(const double rad){

//    const double deg = rad*180/tudat::mathematical_constants::LONG_PI;

    const double deg = rad*180.0/pi;

    return deg;
}


/// B-P frame transformations ///



//! Get transformation quaternion from the Body (B) to the Propulsion (P) frame.
Eigen::Quaterniond getBodyToPropulsionFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd RotationAroundZaxis = Eigen::AngleAxisd(
                -1.0 * thrustAzimuth, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd RotationAroundYaxis = Eigen::AngleAxisd(
                -1.0 * thrustElevation,
                Eigen::Vector3d::UnitY( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( RotationAroundYaxis * RotationAroundZaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}


//! Get transformation matrix from the Body (B) to the Propulsion (P) frame.
Eigen::Matrix3d getBodyToPropulsionFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation )
{
    return getBodyToPropulsionFrameTransformationQuaternion(
            thrustAzimuth, thrustElevation ).toRotationMatrix( );
}

//! Get transformation matrix from the Propulsion (P) to the Body (B) frame.
Eigen::Matrix3d getPropulsionToBodyFrameTransformationMatrix(
    const double thrustAzimuth, const double thrustElevation )
{
    return getBodyToPropulsionFrameTransformationMatrix(
            thrustAzimuth, thrustElevation ).transpose( );
}

//! Get transformation quaternion from the Propulsion (P) to the Body (B) frame.
Eigen::Quaterniond getPropulsionToBodyFrameTransformationQuaternion(
        const double thrustAzimuth, const double thrustElevation )
{
    return getBodyToPropulsionFrameTransformationQuaternion(
            thrustAzimuth, thrustElevation ).inverse( );
}


/// Time functions (from stackoverflow.com/questions/13485266/how-to-have-matlab-tic-toc-in-c and http://en.cppreference.com/w/cpp/chrono/duration/duration_cast)

std::stack<clock_t> tictoc_stack;

typedef std::chrono::high_resolution_clock Clock;
double time1;

auto t1 = Clock::now(); // This is set such that t1 can be changed in the tic() void without it being lost after the void ends

void tic(){ // Creating a similar function to tic toc in Matlab
    tictoc_stack.push(clock());
    t1 = Clock::now();

}


void toc(){
    std::cout<<"CPU time Elapsed: "
            << ((double)(clock() - tictoc_stack.top()))/CLOCKS_PER_SEC<<" seconds"
            << std::endl;
    tictoc_stack.pop();

    auto t2 = Clock::now();
//    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    std::chrono::duration<double> fp_s = t2 - t1;
//    std::cout<<"Time Elapsed more accurate = "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<std::endl;
//    std::cout<<"Time Elapsed even more accurately "<<fp_ms.count()<<" milliseconds"<<std::endl;
    std::cout<<"Wall time Elapsed: "<<fp_s.count()<<" seconds"<<std::endl;
}

Eigen::Vector2d tocNum(){
    double cpuTimeElapsed = ((double)(clock() - tictoc_stack.top()))/CLOCKS_PER_SEC;
    std::cout<<"CPU time Elapsed: "
            << ((double)(clock() - tictoc_stack.top()))/CLOCKS_PER_SEC<<" seconds"
            << std::endl;
    tictoc_stack.pop();

    auto t2 = Clock::now();
//    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    std::chrono::duration<double> fp_s = t2 - t1;
//    std::cout<<"Time Elapsed more accurate = "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<std::endl;
//    std::cout<<"Time Elapsed even more accurately "<<fp_ms.count()<<" milliseconds"<<std::endl;
    std::cout<<"Wall time Elapsed: "<<fp_s.count()<<" seconds"<<std::endl;
    double wallTimeElapsed = fp_s.count();

    Eigen::Vector2d times;
    times(0) = cpuTimeElapsed;
    times(1) = wallTimeElapsed;

    return times;
}
