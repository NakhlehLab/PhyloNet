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

    public String toString () {
	return (map.toString());
    }

    public String toString (NumberFormat numberFormat)
    {
	return (this.toString());
    }
}
