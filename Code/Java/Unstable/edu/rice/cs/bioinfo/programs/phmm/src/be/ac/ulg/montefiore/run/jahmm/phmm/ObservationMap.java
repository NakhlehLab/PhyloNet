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

/**
 * Container for a data column from an MSA.
 * Map from taxon name to observed letter.
 *
 * WARNING - this class is supposed to be immutable.
 * String keys and values are immutable.
 * Use a little extra memory to enforce immutability
 * of underlying Map store.
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

import java.text.NumberFormat;
import java.util.Map;
import java.util.Hashtable;
import be.ac.ulg.montefiore.run.jahmm.Observation;

public class ObservationMap extends Observation
{
    /**
     * The observation's value.
     */
    protected final Map<String,String> map;

    /**
     * An observation that can be described by a string->string map.
     * Although this imposes a memory penalty,
     * force a copy of input map at construction to enforce immutability
     * of this class's objects.
     *
     * @param value The value of this observation.
     */
    public ObservationMap (Map<String,String> inMap)
    {
    map = new Hashtable<String,String>(inMap);
    }

    /**
     * kliu - only expose get method to outside world.
     * Keep this class immutable.
     */
    public String get (String k) {
    return (map.get(k));
    }

    public String toString () {
    return (map.toString());
    }

    public String toString (NumberFormat numberFormat)
    {
    return (this.toString());
    }
}
