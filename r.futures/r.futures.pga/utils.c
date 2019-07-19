/*!
   \file utils.c

   \brief Helper functions

   (C) 2016-2019 by Anna Petrasova, Vaclav Petras and the GRASS Development Team

   This program is free software under the GNU General Public License
   (>=v2).  Read the file COPYING that comes with GRASS for details.

   \author Anna Petrasova
   \author Vaclav Petras
 */

#include <math.h>

/*!
 * \brief Computes euclidean distance in cells (not meters)
 * \param row1 row1
 * \param col1 col1
 * \param row2 row2
 * \param col2 col2
 * \return distance
 */
double get_distance(int row1, int col1, int row2, int col2)
{
    return sqrt((row1 - row2) * (row1 - row2) + (col1 - col2) * (col1 - col2));
}
