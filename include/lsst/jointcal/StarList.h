// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef LSST_JOINTCAL_STAR_LIST_H
#define LSST_JOINTCAL_STAR_LIST_H

#include <string>
#include <list>
#include <iostream>
#include <memory>

#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

class Frame;

//! std::lists of Stars.
/*!     It is a template class, which means that the Star
type remains undefined until a user defines it.
The std::list related operations (insertion,
sort, traversal) are to be carried out using STL
list operations. Most of the
Star operations rely on routines to be provided in
the Star class, usually
user defined. The instanciation of this class for
BaseStar (i.e. the replacement
of the formal parameter 'Star' by 'BaseStar') is
called BaseStarList.
Take care: what is stored is pointers on Star's and
NOT Star's. This implies
that Stars being inserted in the std::list have to be
obtained using 'new'.  */

template <class Star>
class StarList : public std::list<std::shared_ptr<Star>> {
public:
    typedef std::shared_ptr<Star> Element;
    typedef typename std::list<Element>::const_iterator StarCIterator;
    typedef typename std::list<Element>::iterator StarIterator;

    /* constructors */
    //! : default constructor (empty std::list).
    StarList(){};

    /* destructor */
    virtual ~StarList(){};

    //! invokes print(stream) for all Stars in the std::list.
    void print(std::ostream &stream = std::cout) const {
        for (auto &p : *this) p->print(stream);
    }

    //! a model routine to sort the std::list
    /*! see decreasingFlux() to see what it is, if you
       want another sorting criterion) */
    // le premier de la std::liste a le plus grand flux
    void fluxSort();

    //! cuts the end of the std::list
    void cutTail(const int nKeep);

    //! copy the part of the std::list which is included in the frame at the end of another std::list
    void extractInFrame(StarList<Star> &out, const Frame &frame) const;

    //! clears copy and makes a copy of the std::list to copy
    void copyTo(StarList<Star> &copy) const;

    //! Clears the std::list
    void clearList() { cutTail(0); };

    //! enables to apply a geometrical transform if Star is Basestar or derives from it.
    /*! could be extended to other type of transformations. */

    template <class Operator>
    void applyTransform(const Operator &op) {
        for (auto &p : *this) op.transformStar(*(p));
    }
};

//! enables \verbatim  std::cout << my_list; \endverbatim
template <class Star>
std::ostream &operator<<(std::ostream &stream, const StarList<Star> &list) {
    list.print(stream);
    return stream;
}
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_STAR_LIST_H
