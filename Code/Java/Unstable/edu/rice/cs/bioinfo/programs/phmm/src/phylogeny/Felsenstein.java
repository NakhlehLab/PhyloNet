package phylogeny;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

public class Felsenstein {

    static private String[] genes = {"A", "C", "G", "T"};
    static private double lamda = .1;		// this value of the evolutionary constant
	
	 /**
     * On-the-fly version of getLikelihood() implemented for PhyloNet Trees
     * Calculates emission probability 
     * P[O_t | g(s_i), b_{g(s_i)}, \theta ]
     *
     * See writeup for details.
     */
    public double getLikelihoodtree (Tree atree, ObservationMap column) {
	
		double result = 0;
		ArrayList<Double> rootLikelihoods = getLikelihood(atree.getRoot(), column);
		for (int i = 0; i < 4; i++) {
		    // kliu - ArrayList recalculated four times over
		    result += 0.25 * rootLikelihoods.get(i);
		}
	
		// kliu - no need to do this - just overwrite each time
		//clearTree();
	
		return result;
    }
	
    
    /**
     * Jukes-Cantor Model Calculation
     * @param i - A Nucleotide Letter of type char {A, T, G, C}
     * @param j	- A Nucleotide Letter of type char {A, T, G, C}
     * @param u - A Double - base rate substitution/evolutionary constat
     * @param t - A Double - the time interval or branch length time
     * @return returns the Pij value, or the transition probability between two nucleotides given the branch length
     */
    private static double getPij(String i, String j, double u, double t) {
        if (!i.equals(j)) {
            return (0.25 - 0.25 * Math.exp((-4.0 / 3.0) * u * t));
        }
        else
            return (0.25 + 0.75 * Math.exp((-4.0 / 3.0) * u * t));
    }
	
    
    /**
     * @return The node's likelihood arraylist
     */
    public ArrayList<Double> getLikelihood(TNode aNode, ObservationMap column) {
    	ArrayList<Double> likelihood = new ArrayList<Double>();
    	
    	if (aNode.isLeaf()) {
		    //System.out.println(" I am leaf : " + taxa + " and my likelihood array is : " + likelihood);
		    
			String obs = column.get(aNode.getName());
			
			// sets likelihood of leaf based on observation
			int index = -1;
			for (int i = 0; i < genes.length; i++) {
			    if (obs == genes[i]) {
				index = i;
			    }
			}
				
			for (int i = 0; i < 4; i++) {
			    if (i == index) {
				likelihood.add(i, 1.0);
			    }
			    else likelihood.add(i, 0.0);
			}
			
			return likelihood;
		}
		
		
		else {
				
			 Iterator<? extends TNode> childrenIterator= aNode.getChildren().iterator();
			 HashMap<TNode,ArrayList<Double>> childrenLikelihood = new HashMap<TNode,ArrayList<Double>>();
			 ArrayList<TNode> children = new ArrayList<TNode>();
			 
			 while(childrenIterator.hasNext()) {
				 TNode next = childrenIterator.next();
				 children.add(next);
				 childrenLikelihood.put(next, getLikelihood(next, column));
			 }
			 

		    for (int i = 0; i < genes.length; i++) {
				double[] tempvalues = new double[aNode.getChildCount()];
				
				for (int j = 0; j < genes.length; j++) {
					for (int k = 0; k < tempvalues.length; k++) {
						double childLikelihood = childrenLikelihood.get(children.get(k)).get(j);
						tempvalues[k] += childLikelihood * getPij(genes[i], genes[j], lamda, children.get(k).getParentDistance());
					}
				}
				
				double product = 1;
				for (int k = 0; k < tempvalues.length; k++) {
					product *= tempvalues[k];
				}

				likelihood.add(i, product);
		    }
							
		    //System.out.println(" I am an internal node " + this.getTaxa() + " and my likelihood array is : " + likelihood);
	
		    return likelihood;
		}
    }
    
    
}
