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



#ifndef PAGMO_PROBLEM_MAVASCENT_H
#define PAGMO_PROBLEM_MAVASCENT_H

#include <boost/numeric/conversion/cast.hpp>
#include <boost/serialization/export.hpp>
#include <string>

//#include "../tudatBundle/pagmo/src/config.h"
#include "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/config.h"
#include "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/serialization.h"
#include "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/types.h"
#include "/home/stachap/Documents/Thesis/03. Tudat/tudatBundle/pagmo/src/problem/base.h"




namespace pagmo { namespace problem {

/// The MAV ascent problem.
/**
 *
 * This is a box-constrained single-objecive problem.
 * The objective function is the propagation of the ascent trajectory with different possible cost functions;
 * - Minimum mass
 * - Closest correlation to reference ascent trajectory
 *
 */

class __PAGMO_VISIBLE mavAscent : public base
{
        public:
                mavAscent(int = 1);
                base_ptr clone() const;
                std::string get_name() const;
        protected:
                void objfun_impl(fitness_vector &, const decision_vector &) const;
        private:
                friend class boost::serialization::access;
                template <class Archive>
                void serialize(Archive &ar, const unsigned int)
                {
                        ar & boost::serialization::base_object<base>(*this);
                }
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::mavAscent)





#endif // MAVASCENT_H
