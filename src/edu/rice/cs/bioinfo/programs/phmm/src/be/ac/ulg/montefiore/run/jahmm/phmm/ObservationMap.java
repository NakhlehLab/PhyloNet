/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Container for a data column from an MSA.
 * Map from taxon name to observed letter.
 *
 * WARNING - this class is supposed to be immutable.
 * String keys and Character values are immutable.
 * Use a little extra memory to enforce immutability
 * of underlying Map store.
 *
 * Changed from String values to Character values.
 */

package edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.phmm;

import java.text.NumberFormat;
import java.util.Map;
import java.util.Set;
import java.util.Hashtable;
import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm.Observation;

public class ObservationMap extends Observation
{
    /**
     * The observation's value.
     */
    protected final Map<String,Character> map;

    /**
     * An observation that can be described by a string->string map.
     * Although this imposes a memory penalty,
     * force a copy of input map at construction to enforce immutability
     * of this class's objects.
     *
     * @param value The value of this observation.
     */
    public ObservationMap (Map<String,Character> inMap)
    {
	map = new Hashtable<String,Character>(inMap);
    }

    /**
     * kliu - only expose get method to outside world.
     * Keep this class immutable.
     */
    public char get (String k) {
	return (map.get(k).charValue());
    }

    /**
     * Also expose this method.
     */
    public Set<String> keySet () {
	return (map.keySet());
    }

    public String toString () {
	return (map.toString());
    }

    public String toString (NumberFormat numberFormat)
    {
	return (this.toString());
    }
}
