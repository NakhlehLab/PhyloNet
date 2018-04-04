package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.summary;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.Map;

/**
 * Normalized Rooted Branch score (nrBS)
 * Defined in the Appendix of *BEAST paper "Bayesian Inference of Species Trees from Multilocus Data"
 * http://mbe.oxfordjournals.org/content/27/3/570.full
 * Created by wendingqiao on 5/6/16.
 */
public class RFDistanceBL {

    private double _distance;

    public RFDistanceBL(Tree t1, double s1, Tree t2, double s2) {
        _distance = computeDistance(t1, s1, t2, s2);
    }

    public RFDistanceBL(Tree t1, Tree t2) {
        _distance = computeDistance(t1, 1.0, t2, 1.0);
    }

    public double getDistance() {
        return _distance;
    }

    private double computeDistance(Tree t1, double s1, Tree t2, double s2) {
        if(!Trees.leafSetsAgree(t1, t2)) {
            throw new RuntimeException("Trees must have identical leaf sets");
        }
        String[] taxa = t1.getLeaves();
        Map<STITreeCluster, TNode> clusters1 = ((STITree)t1).getClusters(taxa);
        Map<STITreeCluster, TNode> clusters2 = ((STITree)t2).getClusters(taxa);

        double sum = 0;

        for(STITreeCluster cl1: clusters1.keySet()){
            double dist;
            if(clusters2.keySet().contains(cl1)) {
                dist = clusters1.get(cl1).getParentDistance() * s1 - clusters2.get(cl1).getParentDistance() * s2;
            } else {
                dist = clusters1.get(cl1).getParentDistance() * s1;
            }
            sum += dist * dist;
        }
        for(STITreeCluster cl2: clusters2.keySet()){
            if(clusters1.keySet().contains(cl2)) {
                continue; // has been processed
            }
            double dist = clusters2.get(cl2).getParentDistance() * s2;
            sum += dist * dist;
        }
        return Math.sqrt(sum);
    }

    // test
    public static void main(String[] args) {
        String newick1 = "((C:0.002859330148428049,G:0.002859330148428049):0.010473421735219215,(R:0.010936563366414687,(L:0.008340943025951807,(A:0.004491653225077681,Q:0.004491653225077681):0.0038492898008741254):0.00259562034046288):0.0023961885172325784);";
        String newick2 = "((L:0.087137,R:0.087137):0.036282,((G:0.037923,C:0.037923):0.069297,(Q:0.032283,A:0.032283):0.074936):0.016199):0;";
        try{
            NewickReader nr1 = new NewickReader(new StringReader(newick1));
            STITree tree1 = new STITree();
            nr1.readTree(tree1);
            NewickReader nr2 = new NewickReader(new StringReader(newick2));
            STITree tree2 = new STITree();
            nr2.readTree(tree2);
            RFDistanceBL dist = new RFDistanceBL(tree1, 1/0.00286, tree2, 1/0.0379);
            System.out.println(dist.getDistance());
        } catch (Exception ex) {}
    }
}
