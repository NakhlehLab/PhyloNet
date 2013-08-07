/**
 * Miscellaneous utility functions for vectors/matrices.
 * Hmm... might be able to make arbitrarily recursive ArrayList structure.
 * Meh, just use up to 3-d arrays for speed. get() will happen a lot.
 * Later, can parameterize by generic types. Meh, no need for now.
 */

package util;

public class Matrix {
    // base case
    public static boolean checkDimensionsEqual (LogP[] a, LogP[] b) {
	return (a.length == b.length);
    }

    /**
     * Dimensions must match in same row order.
     * Hmm... will this work for multi-dimensional arrays?
     */
    public static boolean checkDimensionsEqual (LogP[][] a, LogP[][] b) {
	if (a.length != b.length) {
	    return (false);
	}

	for (int i = 0; i < a.length; i++) {
	    if (!checkDimensionsEqual(a[i], b[i])) {
		return (false);
	    }
	}

	return (true);
    }

    /**
     * Dimensions must match in same row order.
     * man, really should be recursive version of this for matrices
     * of arbitrary dimension. Requires Java support for 
     * matrices of arbitrary dimension. Couldn't find such support,
     * short of full-blown reflection.
     */
    public static boolean checkDimensionsEqual (LogP[][][] a, LogP[][][] b) {
    	if (a.length != b.length) {
    	    return (false);
    	}

    	for (int i = 0; i < a.length; i++) {
    	    if (!checkDimensionsEqual(a[i], b[i])) {
    		return (false);
    	    }
    	}

	return (true);
    }

    // for set, just use assignment operator to switch references
    // no deep copy needed

    public static String toString (String[] a) {
	String result = "";
	for (String s : a) {
	    result += s + " "; 
	}
	return (result);
    }

    public static String toString (double[] a) {
	String result = "";
	for (double d : a) {
	    result += d + " "; 
	}
	return (result);
    }


    public static String toString (int[] a) {
	String result = "";
	for (int i : a) {
	    result += i + " "; 
	}
	return (result);
    }

    public static String toString (LogP[] a, boolean convertLogFlag) {
	String result = "";
	for (LogP d : a) {
	    result += d.toString(convertLogFlag) + " "; 
	}
	return (result);
    }

    public static String toString (double[][] aa) {
	String result = "";
	for (double[] a : aa) {
	    for (double d : a) { 
		result += d + " "; 
	    }
	    result += "\n";
	}
	return (result);
    }

    public static String toString (int[][] aa) {
	String result = "";
	for (int[] a : aa) {
	    for (int i : a) { 
		result += i + " "; 
	    }
	    result += "\n";
	}
	return (result);
    }

    public static String toString (int[][][] aaa) {
	String result = "";
	int i = 0;
	for (int[][] aa : aaa) {
	    result += i + ": \n"; 
	    result += toString(aa) + "\n";
	    i++;
	}
	return (result);
    }

    public static String toString (LogP[][] aa) {
	String result = "";
	for (LogP[] a : aa) {
	    for (LogP d : a) { 
		result += d + " "; 
	    }
	    result += "\n";
	}
	return (result);
    }

    public static String toString (LogP[][] aa, boolean convertLogFlag) {
	String result = "";
	for (LogP[] a : aa) {
	    result += toString(a, convertLogFlag) + "\n";
	}
	return(result);
    }

    public static String toString (LogP[][][] aaa, boolean convertLogFlag) {
	String result = "";
	int i = 0;
	for (LogP[][] aa : aaa) {
	    result += i + ":" + "\n";
	    result += toString(aa, convertLogFlag) + "\n";
	    i++;
	}
	
	return(result);
    }

    /**
     * Convenience function for dynamically allocating
     * arrays of LogP objects and calling constructor of each
     * array member.
     */
    public static LogP[] createLogPArray (int n) {
	LogP[] a = new LogP[n];
	for (int i = 0; i < a.length; i++) {
	    a[i] = new LogP();
	}
	return (a);
    }

    public static LogP[][] createLogPArray (int n, int o) {
	LogP[][] aa = new LogP[n][o];
	for (int i = 0; i < aa.length; i++) {
	    aa[i] = createLogPArray(o);
	}
	return (aa);
    }

    public static LogP[][][] createLogPArray (int n, int o, int p) {
	LogP[][][] aaa = new LogP[n][o][p];
	for (int i = 0; i < aaa.length; i++) {
	    aaa[i] = createLogPArray(o, p);
	}
	return (aaa);
    }

    public static int sum (int[] ia) {
	int result = 0;
	for (int i : ia) {
	    result += i;
	}
	return (result);
    }

    public static double sum (double[] da) {
	double result = 0.0;
	for (double d : da) {
	    result += d;
	}
	return (result);
    }

    /**
     * convenience normalization functions
     * normalize so entries sum to 1
     */
    public static void normalize (LogP[] p) {
	LogP norm = sum(p);
	
	if (norm.checkNegativeInfinity()) {
	    System.err.println ("ERROR: cannot normalize an array with sum equal to zero. Returning.");
	    return;
	}
	
	for (int i = 0; i < p.length; i++) {
	    p[i].divide(norm);
	}
    }

    public static void normalize (LogP[][] p) {
	LogP norm = sum(p);

	 if (norm.checkNegativeInfinity()) {
	     System.err.println ("ERROR: cannot normalize an array with sum equal to zero. Returning.");
	     return;
	 }

	 for (int i = 0; i < p.length; i++) {
	     for (int j = 0; j < p[i].length; j++) {
		 p[i][j].divide(norm);
	     }
	 }
    }

    
       

    /**
     * return \sum_{i = 0}^{n - 1} a_i
     */
    protected static LogP sum (LogP[] a, boolean[] mask) {
	LogP sum = new LogP(0.0);

	// paranoid
	if (a.length < 0) {
	    System.err.println ("ERROR: empty array in sum(). Returning null.");
	    return (null);
	}

	// paranoid
	if (mask.length != a.length) {
	    System.err.println ("ERROR: mask array length in max() is not equal to input array a length. Returning null.");
	    return (null);
	}

	boolean maskEmptyFlag = true;
	for (boolean b : mask) {
	    if (b) {
		maskEmptyFlag = false;
	    }
	}
	if (maskEmptyFlag) {
	    System.err.println ("ERROR: mask in max() cannot be empty mask! Must select at least one element in input array a. Returning null.");
	    return (null);
	}

	for (int i = 0; i < a.length; i++) {
	    if (mask[i]) {
		sum.add(a[i]);
	    }
	}
	return (sum);
    }

    protected static LogP sum (LogP[] a) {
	LogP sum = new LogP(0.0);

	// paranoid
	if (a.length < 0) {
	    System.err.println ("ERROR: empty array in sum(). Returning null.");
	    return (null);
	}

	for (int i = 0; i < a.length; i++) {
	    sum.add(a[i]);
	}
	return (sum);
    }

    /**
     * Convenience function to return 
     * return \sum_{ij} a_{ij}
     */
    protected static LogP sum (LogP[][] a) {
	LogP sum = new LogP(0.0);

	// paranoid
	if (a.length < 0) {
	    System.err.println ("ERROR: empty array in sum(). Returning null.");
	    return (null);
	}

	for (int i = 0; i < a.length; i++) {
	    sum.add(sum(a[i]));
	}
	return (sum);
    }

    /**
     * no need to do more memory allocation here
     * WARNING! reference returned by this function is to an object in array a!!!
     * return \max_{i has ith position on in bitmask mask} a_i and \argmax of the same expression
     * as [max, argmax]
     *
     * need to define a mask of which array a entries are legal
     * for max check
     */
    protected static void max (LogP[] a, boolean[] mask, Object[] obuffer) {
	// paranoid
	if (a.length < 0) {
	    System.err.println ("ERROR: empty array in max(). Returning.");
	    return;
	}

	// paranoid
	if (mask.length != a.length) {
	    System.err.println ("ERROR: mask array length in max() is not equal to input array a length. Returning.");
	    return;
	}

	boolean maskEmptyFlag = true;
	for (boolean b : mask) {
	    if (b) {
		maskEmptyFlag = false;
	    }
	}
	if (maskEmptyFlag) {
	    System.err.println ("ERROR: mask in max() cannot be empty mask! Must select at least one element in input array a. Returning.");
	    return;
	}

	// paranoid
	if (obuffer.length != 2) {
	    System.err.println ("ERROR: incorrect obuffer dimensions in max(). Returning.");
	    return;
	}

	LogP result = null;
	int resultS = -1;

	for (int i = 0; i < a.length; i++) {
	    if (mask[i] && ((result == null) || a[i].isGreaterThan(result))) {
		result = a[i];
		resultS = i;
	    }
	}

	obuffer[0] = result;
	obuffer[1] = new Integer(resultS);
    }

    /**
     * Convert from arrays of double to arrays of LogP.
     * Assumes that new memory allocated.
     */
    public static void convert (double[] da, LogP[] la) {
	for (int i = 0; i < da.length; i++) {
	    la[i].set(da[i]);
	}
    }

    public static void convert (double[][] daa, LogP[][] laa) {
	for (int i = 0; i < daa.length; i++) {
	    convert(daa[i], laa[i]);
	}
    }
}
