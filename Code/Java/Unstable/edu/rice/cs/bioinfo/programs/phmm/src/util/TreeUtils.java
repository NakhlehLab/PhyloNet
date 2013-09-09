/**
 * Bleh.
 * Isolate name conflict between gsp.ra.Tree and
 * phylonet.....Tree.
 */

package util;

// bleh
//import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
//import gsp.ra.Tree;

public class TreeUtils {
    public static int calculateRobinsonFouldsDistance (edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree t1, 
						       edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree t2) {
	gsp.ra.Tree ut1 = new gsp.ra.Tree();
	gsp.ra.Tree ut2 = new gsp.ra.Tree();
	// kliu - was toNewickWD - I don't think that was correct
	ut1.parseTreeString(t1.toNewick());
	ut2.parseTreeString(t2.toNewick());
	// work with unrooted versions for Robinson-Foulds calculation
	ut1.unroot();
	ut2.unroot();
	return (ut1.computeRobinsonFouldsDistance(ut2));
    }
}
