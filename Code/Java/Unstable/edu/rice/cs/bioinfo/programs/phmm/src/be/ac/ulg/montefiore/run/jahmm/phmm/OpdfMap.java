/* jahmm package - v0.6.1 */

/*
 *  Copyright (c) 2004-2006, Jean-Marc Francois.
 *
 *  This file is part of Jahmm.
 *  Jahmm is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Jahmm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Jahmm; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

import be.ac.ulg.montefiore.run.jahmm.Opdf;
import phylogeny.EvoTree;
import phylogeny.Node;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.HashMap;



/**
 * Utility class for calculating the emission probability
 * of an observation in a PhyloNet-HMM hidden state. 
 * An observation consists of a mapping from taxon names to
 * an observed letter (nucleotide, amino acid, etc.) for a taxon.
 *
 * Under the current PhyloNet-HMM description, emission probabilities
 * consist of two factors: the probability of a gene tree
 * given the parental tree, and the probability of an observed data column
 * given the gene tree.
 * See project writeup for details.
 */
public class OpdfMap
    implements Opdf<ObservationMap>
{
    public static final boolean DEBUG_FLAG = false;

    /**
     * Keep track of containing hidden state.
     * Need to access data associated with the hidden state (gene tree + branch lengths
     * in expected number of substitutions, parental tree + branch lengths in 
     * coalescent units).
     */
    protected HiddenState hiddenState;
    
    /**
     * The base substitution rate.
     */
    protected double lamda = .1;
	
    /**
     * Constructor. For simplicity, pass in reference to the
     * associated hidden state. Need to access data associated with
     * hidden state to perform emission probability calculation on the fly.
     */
    public OpdfMap (HiddenState inHiddenState) {
	setHiddenState(inHiddenState);
    }

    public HiddenState getHiddenState () {
	return (hiddenState);
    }

    public void setHiddenState (HiddenState inHiddenState) {
	this.hiddenState = inHiddenState;
    }
	
    public double probability (ObservationMap o)
    {
	// kliu - on-the-fly emission probability calculation goes here
	// later - add caching dependent on changes to emission probability parameters
	//
	// Emission probability calculation:
	// \begin{eqnarray}
	// e_{s_i, \theta}(O_t) & = & P[O_t | x_t = s_i, \theta] \nonumber \ \
	//  & = & P[O_t | g(s_i), b_{g(s_i)}, \theta ] P[g(s_i), b_{g(s_i)} | T(s_i), c_{T(s_i)}] \nonumber \ \
	//  & = & P[O_t | g(s_i), b_{g(s_i)}, \theta ] P[g(s_i) | T(s_i), c_{T(s_i)}] \nonumber
	// \end{eqnarray}
	//
	// See writeup for details.

	double substitutionModelProbability = hiddenState.getGeneGenealogy().getLikelihood(o);
	// kliu - moved this from emission probability to transition probability
	// see revised writeup
	//double coalescentModelProbability = hiddenState.calculateProbabilityOfGeneGenealogyInParentalTree(DEBUG_FLAG);

	// looks good
	//
	// kliu - testing
	// System.out.println ("Hidden state with parental tree " + hiddenState.getParentalTree().toNewickString() + " and gene genealogy " + hiddenState.getGeneGenealogy().toNewickString() + " :");
	// System.out.println ("substitutionModelProbability: " + substitutionModelProbability);
	// System.out.println ("coalescentModelProbability: " + coalescentModelProbability);

	// kliu - ERROR!
	// This causes emission probability to not be a proper probability distribution!
	// $ \sum_b e_x(b) = \sum_b P[b | g] P[g | T] = P[g | T] \neq 1 $
	// 
	// The problem is the coalescent model contribution. Multiply the non-self-
	// transition probability by the coalescent model factor 
	// and normalize the self-transition probabilities appropriately?
	return (substitutionModelProbability);
    }

    // moved coalescent model contribution to transition probabilities
    //* coalescentModelProbability

    /**
     * TODO. For now, just a warning.
     */
    public ObservationMap generate()
    {	
    	//System.err.println ("ERROR: OpdfMap.generate() not implemented yet. Returning null.");
    	
    	Map<String, String> myMap = new HashMap<String, String>();
    	Node myRoot = hiddenState.geneGenealogy.getRoot();
    	String rootGene= new String();
    	double rand = Math.random();
    	if (rand < 0.25)
    		rootGene = "A";
    	else if (rand < 0.5)
    		rootGene = "T";
    	else if (rand < 0.75)
    		rootGene = "G";
    	else
    		rootGene = "C";
    	
    	// System.out.println("The root gene is: " + rootGene);
    	
    	generateHelper(myMap, myRoot.getLeft(), rootGene);
    	generateHelper(myMap, myRoot.getRight(), rootGene);
    	
    	ObservationMap myColumn = new ObservationMap(myMap);
    	
    	return (myColumn);
    }
    
    
    protected void generateHelper(Map<String, String> yourColumn, Node curNode, String previousGene) {
    	if (curNode.isLeaf()) {
    		String tempTaxa = curNode.getTaxa();
    		String tempGene = selectGene(previousGene, lamda, curNode.getTbranch());
    		yourColumn.put(tempTaxa, tempGene);
    		// System.out.println("Get to the leaf node! " + tempTaxa + ": " + tempGene);
    	}
    	else {
    		String tempGene = selectGene(previousGene, lamda, curNode.getTbranch());
    		// System.out.println("Internal node: " + tempGene);
    		generateHelper(yourColumn, curNode.getLeft(), tempGene);
    		generateHelper(yourColumn, curNode.getRight(), tempGene);
    	}
    }
    
    
    /**
     * Decide which nucleotide a gene will transform to in a time interval of t, according to the
     * base substitution probability.
     * @param geneIn The original gene.
     * @param u The base substitution rate.
     * @param t The time interval
     * @return A randomly selected nucleotide.
     */
    
    protected String selectGene(String geneIn, double u, double t) {
    	String geneOut = new String();
    	ArrayList<String> allGenes = new ArrayList<String>();
    	allGenes.add("A");
    	allGenes.add("T");
    	allGenes.add("G");
    	allGenes.add("C");
    	
    	double rand = Math.random();
    	double p0 = 0.25 + 0.75 * Math.exp((-4.0 / 3.0) * u * t);
    	double p1 = 0.25 - 0.25 * Math.exp((-4.0 / 3.0) * u * t);
    	
    	if (rand < p0) {
    		geneOut = geneIn;
    		return geneOut;
    	}
    	
    	allGenes.remove(geneIn);
    	
    	if (p0 <= rand && rand < p0+p1) {
    		geneOut = allGenes.get(0);
    	}
    	
    	if (p0+p1 <= rand && rand < p0+p1*2) {
    		geneOut = allGenes.get(1);
    	}
    	
    	if (p0+p1*2 <= rand) {
    		geneOut = allGenes.get(2);
    	}
    	
    	return geneOut;
    }
	
	
    public void fit (ObservationMap... oa)
    {
	fit(Arrays.asList(oa));
    }

    public void fit (ObservationMap[] o, double[] weights)
    {
	fit(Arrays.asList(o), weights);
    }

    // kliu - looks like it's just taking the empirical distribution
    // for MLE fit
    //
    // MLE estimate boils down to phylogenetic MLE
    // and coalescent MLE.
    // Here's where we use the hill-climbing, if we so choose (a la RAxML, etc.).
    // Or don't do fitting here - instead iterate between hill-climbing and Baum-Welch E-M.
    //
    // Below re-estimation formula isn't what we want anyways.

    /**
     * MLE estimates for emission probability parameters are obtained
     * via an E-M heuristic incorporating univariate optimization.
     * Estimation occurs outside of Jahmm framework/Baum-Welch E-M.
     *
     * Thus, NOOP here.
     */
    public void fit (Collection<? extends ObservationMap> co)
    {
	// NOOP
	return;
    }
		
    // kliu - looks like just empirical weighted distribution
    // similar to above, but weighted by weights
    public void fit (Collection<? extends ObservationMap> co,
		     double[] weights)
    {
	// NOOP
	return;
    }	
	
    public OpdfMap clone ()
    {	
	return (new OpdfMap(getHiddenState()));
    }
	
	
    public String toString ()
    {
	return toString(NumberFormat.getInstance());
    }
	
	
    public String toString (NumberFormat numberFormat)
    {
	return ("OpdfMap's hiddenState:\n" +
	        hiddenState.toString() + "\n");
    }
	
    protected static final long serialVersionUID = 1L;
}








// public void setProb(int obsNo, double prob)  {
//     this.probabilities[obsNo] = prob;
// }
   
// public void printProbArray() {
//     for (int i = 0; i < probabilities.length; i++)
//         System.out.println(probabilities[i]);
// }


//private static final long serialVersionUID = 1L;
// if (co.isEmpty())
//     throw new IllegalArgumentException("Empty observation set");
		
// for (int i = 0; i < probabilities.length; i++)
//     probabilities[i] = 0.;
		
// for (ObservationInteger o : co)
//     probabilities[o.value]++;
		
// for (int i = 0; i < probabilities.length; i++)
//     probabilities[i] /= co.size();






// if (co.isEmpty() || co.size() != weights.length)
//     throw new IllegalArgumentException();
		
// Arrays.fill(probabilities, 0.);
		
// int i = 0;
// for (ObservationInteger o : co) 
//     probabilities[o.value] += weights[i++];



// if (o.value > probabilities.length-1)
//     throw new IllegalArgumentException("Wrong observation value");
		
// return probabilities[o.value];
