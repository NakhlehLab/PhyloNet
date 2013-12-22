/**
 * Bleh.
 * Isolate name conflict between gsp.ra.Tree and
 * phylonet.....Tree.
 */

package util;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

// bleh
//import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
//import gsp.ra.Tree;

public class TreeUtils {
    /**
     * WARNING - no safety checks.
     * Doesn't check to see if both input trees have the same
     * number of taxa, the same set of taxa, etc.
     */
    public static int calculateRobinsonFouldsDistance (Tree t1, 
						       Tree t2) {
	gsp.ra.Tree ut1 = new gsp.ra.Tree();
	gsp.ra.Tree ut2 = new gsp.ra.Tree();
	// kliu - was toNewickWD - I don't think that was correct
	ut1.parseTreeString(t1.toNewick());
	ut2.parseTreeString(t2.toNewick());

	// work with unrooted versions for Robinson-Foulds calculation
	ut1.unroot();
	ut2.unroot();

	int numInternalEdges1 = ut1.getInternalEdges().length;
	int numInternalEdges2 = ut2.getInternalEdges().length;

	// warning - Robinson-Foulds calculation only works for unrooted
	// trees with at least one internal edge
	// 
	// fails for less than four taxa
	//
	// add guards for these cases
	if ((ut1.getNumTaxa() < 4) || (numInternalEdges1 <= 0)) {
	    return (numInternalEdges2);
	}
	else if ((ut2.getNumTaxa() < 4) || (numInternalEdges2 <= 0)) {
	    return (numInternalEdges1);
	}
	else {
	    return (ut1.computeRobinsonFouldsDistance(ut2));
	}
    }
}