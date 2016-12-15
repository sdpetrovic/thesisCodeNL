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
 *      160921    S.D. Petrovic     File created
 *
 *    References
 *
 *    Notes
 *
 */

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <string>

#include "../tudatBundle/pagmo/src/exceptions.h"
#include "../tudatBundle/pagmo/src/types.h"
#include "../tudatBundle/pagmo/src/problem/base.h"
#include "mavAscent.h"


namespace pagmo { namespace problem {

/// Constructor from dimension.
/**
 * Will construct an n dimensional mavAscent problem.
 *
 * @param[in] n integer dimension of the problem.
 *
 * @see problem::base constructors.
 */
mavAscent::mavAscent(int n):base(n)
{
        // Set bounds.
        set_lb(-500);
        set_ub(500);
}

/// Clone method.
base_ptr mavAscent::clone() const
{
        return base_ptr(new mavAscent(*this));
}

/// Implementation of the objective function.
void mavAscent::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
//        pagmo_assert(f.size() == 1);
//        std::vector<double>::size_type n = x.size();
//        double value=0;

//        for (std::vector<double>::size_type i=0; i<n; i++){
//                value += x[i] * sin(sqrt(fabs(x[i])));
//                }
//                f[0] = 418.9828872724338 * n - value;

//    const decision_vector::size_type n = x.size();
//    f[0]=0;
//    for (decision_vector::size_type i=0; i<n-1; ++i){
//            f[0] += 100 * (x[i]*x[i] -x[i+1])*(x[i]*x[i] -x[i+1]) + (x[i]-1)*(x[i]-1);
//    }
}

std::string mavAscent::get_name() const
{
        return "MAV Ascent";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::mavAscent)
