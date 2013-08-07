/* jahmm package - v0.6.1 */

/*
 *  Copyright (c) 2004-2006, Jean-Marc Francois.
 *
 *  This file is part of Jahmm.
 *  Jahmm is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Jahmm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Jahmm; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */

package be.ac.ulg.montefiore.run.jahmm;

import java.text.*;


/**
 * Observations generated by a Markovian process.  Observations are
 * not mutable.
 */
public abstract class Observation
{
    /**
     * Formats this observation using a number formatter.
     *
     * @param numberFormat A number formatter.
     * @return A string describing this observation.
     */
    public abstract String toString(NumberFormat numberFormat);


    /**
     * Converts this Observation to a {@link java.lang.String String}.
     *
     * @return A descriptive string.
     */
    public String toString()
    {
        return toString(NumberFormat.getInstance());
    }
}

