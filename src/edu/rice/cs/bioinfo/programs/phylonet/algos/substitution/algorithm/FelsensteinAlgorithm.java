package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.JCModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import jeigen.DenseMatrix;
import jeigen.Shortcuts;

import java.io.StringReader;
import java.util.*;

/**
 * This class implements Felsenstein's tree pruning algorithm to computer the likelihood of observations given a gene tree.
 */
public class FelsensteinAlgorithm extends NucleotideProbabilityAlgorithm {
    Tree geneTree;
    SubstitutionModel model;

    /**
     * This is the standard reference for which nucleotide has which index.
     * For example, 'C' has index 1.
     */
    public static final String nucleotides = "ACTG";

    /**
     * Construct an FAlgorithm to run on the passed in data.
     * @param geneTree A gene tree with finite branch lengths.
     * @param model A GTRModel used to calculate the transition probabilities.
     */
    public FelsensteinAlgorithm(Tree geneTree, SubstitutionModel model)
    {
        this.geneTree = geneTree;
        this.model = model;
    }


    /**
     * This creates a column matrix representing the likelihood of a leaf. (IE, the actual nucleiotide at 1, and the rest 0).
     * @param nucleotide A nucleotide from the set 'A','C','T', or 'G'.
     * @return A column matrix
     */
    private static DenseMatrix createLeafMatrix(char nucleotide)
    {
        double [][] temp = {
                {0},
                {0},
                {0},
                {0}
        };

        if (nucleotides.indexOf(nucleotide) == -1)
            throw new RuntimeException(nucleotide +" is not a nucleotide base");
        temp[nucleotides.indexOf(nucleotide)][0] = 1;

        return new DenseMatrix(temp);
    }


    /**
     * This returns the likelihood of a node having each nucleotide
     * @param node The node to find the likelihood of.
     * @param dna A mapping of leaf node names from the gene tree to 'A','C','T', or 'G'.
     * @return A column matrix representing the likelihood. Reference {@link #nucleotides} to see which index corresponds to each base.
     */
    private DenseMatrix probabilityOfNode(TNode node,NucleotideObservation dna)
    {
        if (node.isLeaf())
        {
            return createLeafMatrix(dna.getObservationForAllele(node.getName()));
        }
        else
        {
            DenseMatrix start = Shortcuts.ones(4,1);
            for (TNode child : node.getChildren()) {
                double dist = child.getParentDistance();
                if (Double.isNaN(dist) || Double.isInfinite(dist))
                    throw new IllegalArgumentException("Node : " + child.getName() + " has non-finite dist: "+dist);
                DenseMatrix afterTime = model.getCachedProbabilityMatrix(dist).t();
                start = start.mul(afterTime.mmul(probabilityOfNode(child,dna)));
            }
            return start;
        }
    }



    /**
     * Run the Felstein algorithm and obtain the likelihood of the whole tree.
     * @param dna A mapping of leaf node names from the gene tree to 'A','C','T', or 'G'.
     * @return The likelihood.
     */
    @Override
    public double getProbability(NucleotideObservation dna)
    {
        return probabilityOfNode(geneTree.getRoot(),dna).mul(model.getEquilibriumVector()).sum().s();
    }


    public static double test(String treeString,double val)
    {
        //NewickReader nr = new NewickReader(new StringReader("((c:0.2810896950345153,d:0.29618045148596134):0.8079696326949641,(b:0.5523438948289504,a:0.47318178812839295):0.9983209905780183);"));
        NewickReader nr = new NewickReader(new StringReader(treeString));

        STITree<Double> tree = new STITree<Double>(true);
        try {
            nr.readTree(tree);
        }
        catch(Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        Map<String,String> dnaMap = new HashMap<String,String>();
        dnaMap.put("a","GAATACTGATGATAACGTTTAGTGACAGGTAATCATCGTAGATGGTTCTTTCACTCGTAGAAAAATCGCGATCAATCTCAACAGTTGGGTCAGGTGCAAATGCCGCTCGGGATCATATATAGACCTCTCAGGCTAGAGCCCGATGGGCAAAGATCTGACATCCGCTACCTCCTTTCGCTCCTTTCAGGTGTTCCCTGTCTGATGGATTTGGATCGTGCCACTGCGAATCCCACTGCTCACGCGATGAAATGGTAAGCGTTATACCCGATTTTGCAAAGGCTCGGGTTACATACGGGGGTATGATAGATGCAATACTGGTACTTATGAGACACATTCCAGCCTTTCACCCCTAAAATAGCCGGAGACTCAGCGACGTGGTAGTGTGCTAGTTCTAAAGCAACGAACTGCCTGGGGATGATGAAGAGCGCCTTCCTGTACCCGACTTTAGGCTCGGTACAGCTTCGACACTTTTTAAATCATTGCAATTCAAGGTGATGACATCTGTCCCTAGAAGGAGGCTGGAGGGGTTTCTTAGATTGCACACGTTCCGCAGGCAATAGCTACGCTAGATCTCTTCTTATCGTTACAGTTTGAAGAAATATGAAATCGGATCACGTCTCATTCGTTTGTCAGCTGGCCTGTGGACCTTTATCGTCGCCAAGTCTATTTGGATTGTCACCCCGCAAGTCTGAGCGGGGGACAGTGCATCGCCAGTGGGCACCGTAATCTCATGAGAGGCAGGACTTGTTACGCCTCCGATCTTCAGCCGTGCCGTACGCCGGGGCAAGCGCTTGGGGAATGTGAAATCTATGGCTATGTTTTTGGAACCGGCAGTCGAACTGCATGTAGACGTAGTTTGTATCGACCATGGGAATAGTCGCTAAGGCGACTTTCTTATGCAGACTTGTCGCCAAATATAACACCTCCTACGGCCACGTAAAATCGCCACTTTTAGAGAGAAAGCTGCCTTCTCCCTAGAAGGAAAACTGCAGATGGGTCG");
        dnaMap.put("b","GTATTTTGATGACAACCTTATGGGACAGGGCATCAGCCGGCACGGTTCTTACATTCAGAGCAAAACAACGATCAATCTTAACAAGTTGGTCAGGTTCCAATGCCCCTAGGGGTCGTTTTCAGACCGCGCAACCTAGAGCGTGATGGGCATGAATAGGACCTTCGCTACCGCTTCGCTTTGAATTCAGGTGTTCACTGTTGGTTGAAATTCGATCGTGTTACTGCGGCTCCCGCGGCTTTCCCAACGAGATGTTTAGCGTTATTCCGGACGTTGGAACGGATCGGGTTTCCTACGGGAGTATGACCGAGTCAGAACTAGCACTTATGAGACGCAGTTCGGACTTTCTCCCATACCGAGCCCCGAAAGGGAGAGACGCGGTAGTGTGCTCTATCAAAAGCAACGCACACTATGGCGATGATGCAAAGCGCTCTCGTGTACCCGACCCTTGGAGAGATACAGCTTCGACACTTGTTAAGTCTTTGATGTTCAACGACAGGAGATCTGTCCCTAGAAGGTGGTCGGAGGGCGTTCTTCGATTGCACACTTTCCGCACAACATAGCAACCCCGGATTTCTTCTTTTGCTTACAGTATGCAGACAGACTAAATCTGATCTTGTCTCATTCTTTTGTCACCTGTCCTGAGGAGCTTTCTCGTCCTCCAGTGCATTTGGAATGCGGCGTTCAAATACTGAACCGGGGATACTGTTTCGCTATCGTGCACCTTAACCTCATGAGAGGCAGGACTTTTTACGTCTCCTATCTGCAGCCCTGCCGTACGCAGTGGTGTGCTCAGTGGGACTGTTAAACCTATGATTATGTAGTAAGAACCGTCAGGCGGACGTCGGGTGGATGCAGTTTCTATCAGCTGCAGAAATAGTCTCTTTTGCGACTATTTTATACACACATGGACCCACCTATAACAGCACTTTCGGCCACGTGAGGTCGGATCTTCTAGACGCTTAGCTTCCGTGTCAGTAGAAAGAAAGCTGCTGGTTGATAG");
        dnaMap.put("c","GTAAGCAAATGGCAAAGTCTGCTGATAATTGTGCTGTGGGGTCGCATTGCACTTTCACCCACCCATCCCGAGTACTCCAAAACGGCGGGTTGGTTTCTAAGTGTGCGGGAGTTGATATAGGGATGCGAAGGCCTCGACCGCACAGCCGTAGAATCTCTCACTATCTACCTCCTTCCTCTGCATTGTGGCGTGCCTGGGCTGTCTCATTAAAATCCTCACGACCCAACCAACCCTCGCCAACGTATAAGATGGAGGGCAGTGTCGTGTACGATTGATCGCATTCGGCTTGTACTGCAGGGTTGCTGCAATCGTATTACAGTAGTATTGGATTCAATTTAAGATTTCTGTTCTAAGGAGTACCGCGATTGAGCGGCACGGCCGTTTCCAGGGTCTTAAAAAGCACAGATCACTGCCCTTATTCACAGAATCCTCCTTGTCGCTGCCCAAGTATGACTAAAGCTATGATCCTTGCGTTATTGTTGAGAGTTGCCGCCATGAGATAAATAATTTGACGGGGGTCGCTGGTGCGTATTACAGTAGACATATTCCAGATACAATAGCCAAATTGCAATTCTTGTAAGGCTACTGGCCTCAAGAGATATGCGATCGGAGGGTGTACCCCCAATTTCTGTAACGTTCTGGGACTATTTGCCTGCGCCCCGTATATGGTGAGGTACGTCCTGCAATTTGTAAGGGAGGAGTGTTGAGCTTGCCTGAGCCCCTTAACTGACATGAAAGAGAGAATTATTACTACTCCTATCGTTCGCCTACTCGAAGAAATGGGCCTTCTGCAGGCCACTGGAATGCGGAGGTCTACAGGTGATGCAGGTTACATAGGACGTGGTTGTGAGGAGGCTGATCTCCACCGGCTAAATAGTTACTTGGGAGGGTTCTTCATACACGATTGTTGATGAGTATATCAATGACTAGATAATCATGAGGACGCCTGTTTTAGGCCGTAGGTAATGTTTTCAGCACTGCGTAGACTAGACACCGGTGG");
        dnaMap.put("d","CTAAGCAACTGGCAAAGTCTGCTGATGCTTGCGCGGTGGCGGCGCGTAGCTCTTTCACCCAACCCTGTCGAGGACTCCAAAACCGCGGGGTGGTTTCTAAATGTGCCGGCGCTAATGCATGGATCCCAAGGCCTAGAGCGCACAGCCGTAGAATTTCTCCCTAGCTACCTCCTTCCTCGGCATTGGGGCGTGACAGGGCCGTCCTTCTAAGATCGTGACGACCCAACCAACCCTCGCCACCATATCCGATGGAGGGCAGCGTCGTGTACGTTTGACCGCACTCGGCTTGTGCTACAGGTTTGCTGCAATCGCATTACGGCGTTATTAGATTCAATTGAAGCTTTCTATTGTAACGGATACCGCGATTGGTCGGCGCGGCTGTTTCGAGGGTATTGAAAAGAACAGATCACGGTCCTTATTCATAGATTCCTCGTTGTCGCTACCCTAGAATGACTAAAGCTATGATCCTCCCGTTATTGTTGACCATGGACGTCATGGGATATATTATTTAACGGGGGTCGCTGGTGCTTAATACATTAGATATATTCCAGATACAATGGCCAAGTAGCAATGCTTCTAAAGCTGCTGGCCTCGAGAGATACGCAATGGGTTGATGGAGCCCTTTTTTCTCTACCGCTTCGGGAGGATTTGCCGGCGCCCCGTAGATGGTGAGGGACGTCCTGCGATTTGAAAGGGAGTAGTGTTGAGCTTGACCGCGCCCCATAAGTGAAACGAAAGAGAGAAGTATTAAGACTCCTATCTTTGGCCCACTCGAAGAAATGGGCCTTCTGCAGGCCACTGCAAAGCGTAGTTCTACGTGTGAGGGAGGTTCAATAGTACATGGTTGTGAGGAAGTTGATCTGAACCGGATACATAGTGAATTGCGAGCGTTCTTCATACACTATTGTTGCAGATTATATCAATGCCTAGATAATCATGAAGAGGCCTGGTTTAGGCCGTAGGTAAAGTTTTCAGCACCTCGAAGACTATATTCCGTTGG");


        double[] e = {.1,.2,.3,.4};
        double[] t = {val,.9,.3,.2,.5,.2};
        SubstitutionModel model = new JCModel(e,t);

        double sum = 0;
        for (NucleotideObservation nucleotideObservation : NucleotideObservation.getObservations(dnaMap)) {
            //System.out.println(nucleotideObservation.toString());
            sum += Math.log(new FelsensteinAlgorithm(tree,model).getProbability(nucleotideObservation));
        }
        return sum;
    }


    public static void main(String[] args)
    {
        //System.out.println("Before:");
        //System.out.println(test("((c:0.2810896950345153,d:0.29618045148596134):0.8079696326949641,(b:0.5523438948289504,a:0.47318178812839295):0.9983209905780183);",.4889));
        //System.out.println(test("((c:0.1,d:0.1):0.4,(b:0.2,a:0.2):.3);",1));
        //System.out.println("After:");

        NewickReader nr = new NewickReader(new StringReader("((SA:1.0,SB:1.0)SAB:0.5,(SC:1.0,SD:1.0)SCD:0.5)root"));

        STITree<Double> tree = new STITree<Double>(true);
        try {
            nr.readTree(tree);
        }
        catch(Exception e) {
            e.printStackTrace();
            return;
        }

        Map<String,Character> omap = new HashMap<String,Character>();
        omap.put("SA", 'G');
        omap.put("SB", 'C');
        omap.put("SC", 'T');
        omap.put("SD", 'G');

        double[] e = {.1,.2,.3,.4};
        double[] t = {.4,.9,.3,.2,.5,.2};
        SubstitutionModel model = new GTRModel(e,t);

        System.out.println("Debug:");
        System.out.println("One:");
        System.out.println(new FelsensteinAlgorithm(tree,model).getProbability(new OneNucleotideObservation(omap)));
        System.out.println("Two");
        System.out.println(new Felsenstein().getLikelihoodtree(tree,omap,model));


    }
}



