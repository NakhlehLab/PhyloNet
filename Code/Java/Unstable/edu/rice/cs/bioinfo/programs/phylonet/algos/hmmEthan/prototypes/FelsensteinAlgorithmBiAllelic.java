package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.prototypes;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.NucleotideProbabilityAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import jeigen.DenseMatrix;
import jeigen.Shortcuts;

import java.io.StringReader;
import java.util.HashMap;
import java.util.Map;

/**
 * This class implements Felstines tree pruning algorithm to computer the likelihood of observations given a gene tree.
 */
public class FelsensteinAlgorithmBiAllelic extends NucleotideProbabilityAlgorithm {
    Tree geneTree;
    SubstitutionModel model;

    /**
     * This is the standard reference for which nucleotide has which index.
     * For example, 'C' has index 1.
     */
    static String nucleotides = "01";

    /**
     * Construct an FAlgorithm to run on the passed in data.
     *
     * @param geneTree A gene tree with finite branch lengths.
     * @param model    A GTRModel used to calculate the transition probabilities.
     */
    public FelsensteinAlgorithmBiAllelic(Tree geneTree, SubstitutionModel model) {
        this.geneTree = geneTree;
        this.model = model;
    }


    /**
     * This creates a column matrix representing the likelihood of a leaf. (IE, the actual nucleiotide at 1, and the rest 0).
     *
     * @param nucleotide A nucleotide from the set 'A','C','T', or 'G'.
     * @return A column matrix
     */
    private static DenseMatrix createLeafMatrix(char nucleotide) {
        double[][] temp = {
                {0},
                {0},
        };

        if (nucleotides.indexOf(nucleotide) == -1)
            throw new RuntimeException(nucleotide + " is not a BiAllelic base");
        temp[nucleotides.indexOf(nucleotide)][0] = 1;

        return new DenseMatrix(temp);
    }

    /**
     * This returns the likelihood of a node having each nucleotide
     *
     * @param node The node to find the likelihood of.
     * @param dna  A mapping of leaf node names from the gene tree to 'A','C','T', or 'G'.
     * @return A column matrix representing the likelihood. Reference {@link #nucleotides} to see which index corresponds to each base.
     */
    private DenseMatrix probabilityOfNode(TNode node, NucleotideObservation dna) {
        if (node.isLeaf()) {
            return createLeafMatrix(dna.getObservationForAllele(node.getName()));
        } else {
            DenseMatrix start = Shortcuts.ones(2, 1);
            for (TNode child : node.getChildren()) {
                double dist = child.getParentDistance();
                if (Double.isNaN(dist) || Double.isInfinite(dist))
                    throw new IllegalArgumentException("Node : " + child.getName() + " has non-finite dist: " + dist);
                DenseMatrix afterTime = model.getProbabilityMatrix(dist).t();
                start = start.mul(afterTime.mmul(probabilityOfNode(child, dna)));
            }
            return start;
        }
    }


    Map<NucleotideObservation, Double> memo = new HashMap<NucleotideObservation, Double>();

    /**
     * Run the Felstein algorithm and obtain the likelihood of the whole tree.
     *
     * @param dna A mapping of leaf node names from the gene tree to 'A','C','T', or 'G'.
     * @return The likelihood.
     */
    @Override
    public double getProbability(NucleotideObservation dna) {
        Double result = memo.get(dna);
        double actualResult;

        if (result == null) {
            actualResult = probabilityOfNode(geneTree.getRoot(), dna).mul(model.getEquilibriumVector()).sum().s();
            memo.put(dna, actualResult);
            return actualResult;

        } else
            return result;
    }




    public static void main(String[] args) {
        NewickReader nr = new NewickReader(new StringReader("((SA:1.0,SB:1.0)SAB:0.5,(SC:1.0,SD:1.0)SCD:0.5)root"));

        STITree<Double> tree = new STITree<Double>(true);
        try {
            nr.readTree(tree);
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        Map<String, Character> omap = new HashMap<String, Character>();
        omap.put("SA", 'G');
        omap.put("SB", 'C');
        omap.put("SC", 'T');
        omap.put("SD", 'G');

        double[] e = {.1, .2, .3, .4};
        double[] t = {.4, .9, .3, .2, .5, .2};
        SubstitutionModel model = new GTRModel(e, t);

        System.out.println(new FelsensteinAlgorithmBiAllelic(tree, model).getProbability(new OneNucleotideObservation(omap)));
    }
}

