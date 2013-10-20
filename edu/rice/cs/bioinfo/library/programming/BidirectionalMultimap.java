/**
 * BidirectionalMultimap.java
 *
 * Maintains many-to-many relationship between key set and value set.
 * Don't bother with interface, etc. For now, just assume two
 * backing Hashtables, with all of their constraints.
 *
 * WARNING! Both K and V *MUST* have meaningful equals() and hashCode()
 * implementations! (Similarly to Hashtable specification.)
 */

package edu.rice.cs.bioinfo.library.programming;

import java.util.Map;
import java.util.Set;
import java.util.Hashtable;
import java.util.HashSet;

public class BidirectionalMultimap<K,V> 
    // Map not applicable here, since
    // both directions are multimaps.
    //implements Map<K,V> 
    // Also, Dictionary<K,V> abstract class is obsolete.
{
    // K' \subset type K
    // V' \subset type V
    // Multimap - key k from K' -> set of values from V'
    protected Hashtable<K,HashSet<V>> map;
    // Reversed multimap - value v from V' -> set of keys k from K'
    // s.t. for all k, map.get(k) includes v
    protected Hashtable<V,HashSet<K>> rmap;

    /**
     * Create an empty BidirectionalMultimap.
     */
    public BidirectionalMultimap () {
	// Don't bother with interface yet.
	// Assume hash-based backing stores.
	map = new Hashtable<K,HashSet<V>>();
	rmap = new Hashtable<V,HashSet<K>>();
    }

    // Don't know why type system is complaining.
    // Too much of a PITA. Change to Hashtable and HashSet everywhere.
    // Probably related to type erasure implementation of generic types.
    /**
     * Isolate implementation-specific details as much as possible.
     */
    protected HashSet<K> createNewKeySet () {
	return (new HashSet<K>());
    }

    protected HashSet<V> createNewValueSet () {
	return (new HashSet<V>());
    }

    /**
     * Number of keys k in K'.
     */
    public int sizeKeys () {
	return (map.size());
    }

    /**
     * Number of values v in V'.
     */
    public int sizeValues() {
	return (rmap.size());
    }

    /**
     * Returns set of all keys.
     */
    public Set<K> keys () {
	return (map.keySet());
    }
    
    /**
     * Returns set of all values.
     */
    public Set<V> values () {
	return (rmap.keySet());
    }

    public boolean containsKey (Object key) {
	return (map.containsKey(key));
    }

    public boolean containsValue (Object value) {
	return (rmap.containsKey(value));
    }

    /**
     * Checks if number of entries is zero.
     */
    public boolean isEmpty () {
	// paranoid
	if (map.isEmpty() ^ rmap.isEmpty()) {
	    throw (new RuntimeException("ERROR: one of map/rmap is empty while the other isn't."));
	}
	return (map.isEmpty());
    }

    /**
     * Adds a (k,v) \in (K',V') pair into
     * the BidirectionalMultimap.
     *
     * Since get() doesn't follow Map contract,
     * put() also doesn't - no return value for put().
     */
    public void put (K key, V value) {
	// add to forward multimap
	if (!map.containsKey(key)) {
	    map.put(key, createNewValueSet());
	}
	HashSet<V> sv = map.get(key);
	sv.add(value);

	// add to reverse multimap
	if (!rmap.containsKey(value)) {
	    rmap.put(value, createNewKeySet());
	}
	HashSet<K> sk = rmap.get(value);
	sk.add(key);
    }

    /**
     * Retrieve the set of values mapping to a key k.
     */
    public Set<V> get (K key) {
	return (map.get(key));
    }

    /**
     * Retrieve the set of keys mapping to a value v.
     * Won't compile vs. above in Java due to generics support through type erasure.
     */
    public Set<K> rget (V value) {
	return (rmap.get(value));
    }

    /**
     * Is the relation a function?
     */
    public boolean checkFunctional () {
	for (K key : keys()) {
	    if (get(key).size() > 1) {
		return (false);
	    }
	}

	return (true);
    }

    /**
     * Does the relation have the injective property?
     */
    public boolean checkInjective () {
	for (V value : values()) {
	    if (rget(value).size() > 1) {
		return (false);
	    }
	}

	return (true);
    }

    /**
     * Is the relation a bijective function, 
     * i.e. does the relation have both the functional and injective properties 
     * (since surjective property is implied)?
     */
    public boolean checkBijective () {
	return (checkFunctional() && checkInjective());
    }

    /**
     * No support for remove() operation yet.
     */

    /**
     * For diagnostic purposes.
     */
    public String toString () {
	String s = "";
	s+= "Keys:\n";
	for (K k : keys()) {
	    s += k.toString() + "\n";
	}

	s+= "Values:\n";
	for (V v : values()) {
	    s += v.toString() + "\n";
	}

    	s += "Forward:\n";
    	for (K k : keys()) {
    	    for (V v : get(k)) {
    		s += k.toString() + ": " + v.toString() + "\n";
    	    }
    	}

    	s += "Reverse:\n";
    	for (V v : values()) {
    	    for (K k : rget(v)) {
    		s += k.toString() + ": " + v.toString() + "\n";
    	    }
    	}

	return (s);
    }

    // Looks good.
    // public static void test () {
    // 	Hashtable<LengthParameter,HashSet<Integer>> ht = new Hashtable<LengthParameter,HashSet<Integer>>();
    // 	LengthParameter t = new LengthParameter("t0",3.0);
    // 	ht.put(t, new HashSet<Integer>());
    // 	for (Map.Entry entry : ht.entrySet()) {
    // 	    System.out.println ("key: |" + entry.getKey() + "| value: |" + entry.getValue() + "| contains: |" + ht.containsKey(entry.getKey()) + "|");
    // 	}
    // 	// this also doesn't work???
    // 	System.out.println ("where? " + t + " " + ht.size() + " " + ht.containsKey(t));



    // 	// Hashtable<String,HashSet<Integer>> h = new Hashtable<String,HashSet<Integer>>();
    // 	// h.put("foobar", new HashSet<Integer>());
    // 	// System.out.println ("h: " + h.containsKey("foobar"));



    // 	BidirectionalMultimap<LengthParameter,Integer> bm = new BidirectionalMultimap<LengthParameter,Integer>();
    // 	bm.put(new LengthParameter("t1",3.0), new Integer(0));
    // 	bm.put(new LengthParameter("t2",4.0), new Integer(1));
    // 	bm.put(new LengthParameter("t3",5.0), new Integer(2));

    // 	System.out.println ("Forward:");
    // 	for (LengthParameter k : bm.keys()) {
    // 	    for (Integer v : bm.get(k)) {
    // 		System.out.println (k + ": " + v);
    // 	    }
    // 	}

    // 	System.out.println ("Reverse:");
    // 	for (Integer v : bm.values()) {
    // 	    for (LengthParameter k : bm.rget(v)) {
    // 		System.out.println (k + ": " + v);
    // 	    }
    // 	}

    // 	System.out.println ("Functional? " + bm.checkFunctional());
    // 	System.out.println ("Injective? " + bm.checkInjective());
    // 	System.out.println ("Bijective? " + bm.checkBijective());
    // }

    // public static void main (String[] args) {
    // 	test();
    // }








   


    // /**
    //  * Just need some storage for constraint.
    //  * Immutable after creation.
    //  */
    // public static class LengthParameterConstraintSet {
    // 	protected String name;
    // 	protected double value;

    // 	public LengthParameterConstraintSet (String inName, double inValue) {
    // 	    setName(inName);
    // 	    setValue(inValue);
    // 	}

    // 	public String getName () {
    // 	    return (name);
    // 	}

    // 	public double getValue () {
    // 	    return (value);
    // 	}

    // 	// immutable once created
    // 	protected void setName (String inName) {
    // 	    name = inName;
    // 	}

    // 	protected void setValue (double inValue) {
    // 	    value = inValue;
    // 	}

    // 	// LengthParameterConstraintSet objects are hashCode() and equals() equivalent to
    // 	// their name Strings
	
    // 	public String toString () {
    // 	    return (name.toString());
    // 	}

    // 	public int hashCode () {
    // 	    return (name.hashCode());
    // 	}

    // 	// derp
    // 	public boolean equals (Object obj) {
    // 	    if (!(obj instanceof LengthParameterConstraintSet)) {
    // 		return (false);
    // 	    }
    // 	    LengthParameterConstraintSet lpcs = (LengthParameterConstraintSet) obj;
    // 	    return (name.equals(lpcs.getName()));
    // 	}
    // }
    
    // /**
    //  * Simple state container for length parameter.
    //  * Reduce indirection between IDs and objects a bit.
    //  *
    //  * Don't clone.
    //  */
    // public static class LengthParameter extends LengthParameterConstraintSet {
    // 	public LengthParameter (String inName, double inValue) {
    // 	    super(inName, inValue);
    // 	}

    // 	public void setName (String inName) {
    // 	    super.setName(inName);
    // 	}

    // 	public void setValue (double inValue) {
    // 	    super.setValue(inValue);
    // 	}

    // 	// derp
    // 	public boolean equals (Object obj) {
    // 	    if (!(obj instanceof LengthParameter)) {
    // 		return (false);
    // 	    }
    // 	    LengthParameter lp = (LengthParameter) obj;
    // 	    return (name.equals(lp.getName()));
    // 	}
    // }


}
