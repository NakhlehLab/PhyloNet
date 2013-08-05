package edu.rice.cs.bioinfo.library.programming;

import java.util.Map;
import java.util.Set;
import java.util.Hashtable;
import java.util.HashSet;

/**
 * Maintains a bijective map.
 * Disallow null key/null value.
 */

public class BijectiveHashtable<K,V> {
    protected Hashtable<K,V> map;
    protected Hashtable<V,K> rmap;
    
    public BijectiveHashtable () {
	map = new Hashtable<K,V>();
	rmap = new Hashtable<V,K>();
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

    public void put (K key, V value) {
	// collisions need to obliterate entries in both maps
	// max two deletions
	// warn in this case
	// may not be desired behavior
	// deletions in case
	// of collision for either key or value
	if (containsKey(key)) {
	    System.err.println ("WARNING: BijectiveHashtable overwriting existing entry " + key + "->" + get(key) +".");
	    rmap.remove(get(key));
	    map.remove(key);
	}

	if (containsValue(value)) {
	    System.err.println ("WARNING: BijectiveHashtable overwriting existing entry " + rget(value) + "->" + value +".");
	    map.remove(rget(value));
	    rmap.remove(value);
	}

	map.put(key, value);
	rmap.put(value, key);
    }

    public V get (K key) {
	return (map.get(key));
    }

    /**
     * la raison d etre
     */
    public K rget (V value) {
	return (rmap.get(value));
    }

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
	    s += k.toString() + ": " + get(k).toString() + "\n";
    	}

    	s += "Reverse:\n";
    	for (V v : values()) {
	    s += rget(v).toString() + ": " + v.toString() + "\n";
    	}

	return (s);
    }

    // Looks good.
    public static void test () {
    	BijectiveHashtable<String,Integer> bm = new BijectiveHashtable<String,Integer>();
    	bm.put("t1", new Integer(0));
    	bm.put("t2", new Integer(1));
    	bm.put("t3", new Integer(2));

	System.out.println (bm);

    	bm.put("t3", new Integer(42));
    	bm.put("foobar", new Integer(0));
    	bm.put("t5", new Integer(101));

	System.out.println (bm);
    }

    public static void main (String[] args) {
    	test();
    }


}
