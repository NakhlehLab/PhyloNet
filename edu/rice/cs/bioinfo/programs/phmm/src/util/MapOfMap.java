/**
 * Silly class. Only reason this exists is because
 * we want to automatically create empty instances of 
 * the inner map/check for non-nullness of inner map.
 *
 * Use a hashtable for now.
 */

package util;

import java.util.Map;
import java.util.Set;
import java.util.Hashtable;

public class MapOfMap<K1,K2,V> {
    Map<K1,Map<K2,V>> mmap;

    public MapOfMap () {
	mmap = new Hashtable<K1,Map<K2,V>>();
    }

    public boolean isEmpty () {
	return (mmap.isEmpty());
    }

    public void clear () {
	for (K1 k1 : mmap.keySet()) {
	    if (mmap.get(k1) != null) {
		mmap.get(k1).clear();
	    }
	}
	mmap.clear();
    }

    public void clear (K1 k1) {
	if (mmap.containsKey(k1)) {
	    if (mmap.get(k1) != null) {
		mmap.get(k1).clear();
	    }
	    mmap.remove(k1);
	}
    }

    public boolean contains (K1 k1, K2 k2) {
	return (mmap.containsKey(k1) &&
		(mmap.get(k1) != null) &&
		mmap.get(k1).containsKey(k2));
    }

    public void put (K1 k1, K2 k2, V v) {
	if (!mmap.containsKey(k1)) {
	    mmap.put(k1, new Hashtable<K2,V>());
	}
	
	mmap.get(k1).put(k2, v);
    }

    public V get (K1 k1, K2 k2) {
    	if (!contains(k1, k2)) {
    	    System.err.println ("ERROR: keys " + k1 + " and " + k2 + " are unknown in MapOfMap object. Returning null to signal error.");
	    return (null);
    	}
	
    	return (mmap.get(k1).get(k2));
    }

    /**
     * Returns all values for first parameter of MapOfMap
     */ 
    public Set<K1> firstKeySet () {
	return (mmap.keySet());
    }
    
    public Map<K2,V> get (K1 k1) {
	// returns null if k1 unknown in MapOfMap object
	return (mmap.get(k1));
    }

    public String toString () {
	StringBuffer sb = new StringBuffer();
	for (K1 k1 : mmap.keySet()) {
	    if (mmap.get(k1) == null) {
		System.err.println ("ERROR: value " + k1 + " for first argument in MapOfMap maps to null entry. Returning empty string to signal error.");
		return ("");
	    }

	    for (K2 k2 : mmap.get(k1).keySet()) {
		sb.append(k1.toString() + " -> " + k2.toString() + " -> " + this.get(k1, k2).toString() + "\n");
	    }
	}
	
	return (sb.toString());
    }
}
