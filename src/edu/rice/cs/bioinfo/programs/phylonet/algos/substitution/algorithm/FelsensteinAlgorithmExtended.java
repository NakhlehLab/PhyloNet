package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm;


import java.io.StringReader;
import java.util.HashMap;
import java.util.Map;

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


public class FelsensteinAlgorithmExtended extends NucleotideProbabilityAlgorithm {
    Tree geneTree;
    SubstitutionModel model;
	    

	public static final String nucleotides = "ACTG";
	
	public FelsensteinAlgorithmExtended(Tree geneTree, SubstitutionModel model)
    {
        this.geneTree = geneTree;
        this.model = model;
    }
	
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

	/*
	 * Obtain the likelihood of the whole tree.
	 * */
    public double getProbability(NucleotideObservation dna)
    {
        return probabilityOfNode(geneTree.getRoot(),dna).mul(model.getEquilibriumVector()).sum().s();
    }

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
                DenseMatrix afterTime = model.getProbabilityMatrixIntegrated().t();
                start = start.mul(afterTime.mmul(probabilityOfNode(child,dna)));
            }
            
            return start;
        }
    }
	
	
    /*
     *  Test case for Jukes-Cantor.
     */

    public static void runJCtest() {
    	//Original Felsenstein's need branch lengths.		
		NewickReader nr = new NewickReader(new StringReader("((SA: 1.0,SB: 1.0)SAB: 1.0, SC: 1.0)root"));
		STITree<Double> tree = new STITree<Double>(true);
        try {
            nr.readTree(tree);
        }
        catch(Exception e) {
            e.printStackTrace();
            return;
        }
 
        // Create the Jukes-Cantor model with mu = 0.1.
        double[] e = {.1,.2,.3,.4};
        double[] t = {0.1,.9,.3,.2,.5,.2};
        SubstitutionModel model = new JCModel(e,t);
        
        
        
        Map<String,Character> omap = new HashMap<String,Character>();

        
        double sumResults = 0.0;
        for (int i = 0; i < nucleotides.length(); i++) {
        	for (int j = 0; j < nucleotides.length(); j++) {
        		for(int k = 0; k < nucleotides.length(); k++) {
        	        omap.put("SA", nucleotides.charAt(i));
        	        omap.put("SB", nucleotides.charAt(j));
        	        omap.put("SC", nucleotides.charAt(k));
        	        
        	        double result =  new FelsensteinAlgorithmExtended(tree,model).getProbability(new OneNucleotideObservation(omap));
        	        sumResults += result;
        	        
        	        System.out.println("Pattern:" + nucleotides.charAt(i)+nucleotides.charAt(j)+nucleotides.charAt(k) + "  Probability:" + result);
        		}
        	}
        }
        
        
        // Maximum allowed error in calculations involving 'double'.
        // Ideally all of the results should sum to 1.
        double epsilon = 0.0000001;
        if (Math.abs(sumResults - 1.0) > epsilon) {
        	System.out.println("ERROR: Error in precision!");
        }
        else {
        	System.out.println("Success: Probabilities of all possible trees sum to 1.0 considering precision.");
        	System.out.println("Actual sum of all probabilities = " + sumResults);
        }   	
    }
    
    
    /*
     *  Test case for GTR.
     */

    public static void runGTRtest() {
    	//Original Felsenstein's need branch lengths.		
		NewickReader nr = new NewickReader(new StringReader("((SA: 1.0,SB: 1.0)SAB: 1.0, SC: 1.0)root"));
		STITree<Double> tree = new STITree<Double>(true);
        try {
            nr.readTree(tree);
        }
        catch(Exception e) {
            e.printStackTrace();
            return;
        }
 
        // Create the GTR model with 4+6 parameters.
        double[] e = {.1,.2,.3,.4};
        double[] t = {0.1,.9,.3,.2,.5,.2};
        SubstitutionModel model = new GTRModel(e,t);
        

        Map<String,Character> omap = new HashMap<String,Character>();
  
        double sumResults = 0.0;
        for (int i = 0; i < nucleotides.length(); i++) {
        	for (int j = 0; j < nucleotides.length(); j++) {
        		for(int k = 0; k < nucleotides.length(); k++) {
        	        omap.put("SA", nucleotides.charAt(i));
        	        omap.put("SB", nucleotides.charAt(j));
        	        omap.put("SC", nucleotides.charAt(k));
        	        
        	        double result =  new FelsensteinAlgorithmExtended(tree,model).getProbability(new OneNucleotideObservation(omap));
        	        sumResults += result;
        	        
        	        System.out.println("Pattern:" + nucleotides.charAt(i)+nucleotides.charAt(j)+nucleotides.charAt(k) + "  Probability:" + result);
        		}
        	}
        }
        
        
        // Maximum allowed error in calculations involving 'double'.
        // Ideally all of the results should sum to 1.
        double epsilon = 0.0000001;
        if (Math.abs(sumResults - 1.0) > epsilon) {
        	System.out.println("ERROR: Error in precision!");
        }
        else {
        	System.out.println("Success: Probabilities of all possible trees sum to 1.0 considering precision.");
        	System.out.println("Actual sum of all probabilities = " + sumResults);
        }   	
    }
        
    
    
    /*
     *  Main runs couple of test cases that must always pass,
     *  irrespective of the above algorithm.
     */
    public static void main(String[] args)
    {
    	runJCtest();
    	runGTRtest();
	}
	
}