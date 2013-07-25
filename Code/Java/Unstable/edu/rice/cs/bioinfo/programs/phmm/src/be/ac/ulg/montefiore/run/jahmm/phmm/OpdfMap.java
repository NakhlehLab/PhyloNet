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
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.Vector;
import java.util.List;
import java.io.StringReader;
import java.io.IOException;

// kliu - Phylonet support libraries
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;



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
    /**
     * Keep track of containing hidden state.
     * Need to access data associated with the hidden state (gene tree + branch lengths
     * in expected number of substitutions, parental tree + branch lengths in 
     * coalescent units).
     */
    protected HiddenState hiddenState;

    /**
     * For coalescent model calculations.
     */
    protected GeneTreeProbability gtp;
	
    /**
     * Constructor. For simplicity, pass in reference to the
     * associated hidden state. Need to access data associated with
     * hidden state to perform emission probability calculation on the fly.
     */
    public OpdfMap (HiddenState inHiddenState) {
	setHiddenState(inHiddenState);
	gtp = new GeneTreeProbability();
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
	// TODO - add in coalescent model contribution
	//
	// later - add caching dependent on changes to emission probability parameters
	return (hiddenState.getGeneGenealogy().getLikelihood(o));
    }

    /**
     * Perform standard coalescent model calculation to obtain
     * probability P[g(s_i) | T(s_i), c_{T(s_i)}] of observing a gene genealogy given a parental tree.
     *
     * See writeup for details.
     *
     * Might move this later to HiddenState.
     *
     * Use code from ComputeGTProb.
     * WARNING - returns likelihood, *NOT* log likelihood!
     */
    protected double calculateProbabilityOfGeneGenealogyInParentalTree (boolean debugFlag) {
	Network<Double> parentalTree = convertPHMMTreeToPhyloNetNetwork(hiddenState.getParentalTree());
	Map<String,String> alleleToSpeciesMapping = hiddenState.getAlleleToSpeciesMapping();
	Tree geneGenealogy = convertPHMMTreeToPhyloNetTree(hiddenState.getGeneGenealogy());
	// A list with one element. Inefficient - consider doing a one-shot approach later.
	Vector<Tree> geneGenealogies = new Vector<Tree>();
	geneGenealogies.add(geneGenealogy);

	gtp.emptyState();

	// calculation under model from Yu et al. 2012
	// crap - this method requires Network<Double>
	// Yun uses Double to store hybridization probabilities during calculation
        List<Double> probList = gtp.calculateGTDistribution(parentalTree, geneGenealogies, alleleToSpeciesMapping, debugFlag);

	// should only be a single entry
	if (probList.size() != 1) {
	    System.err.println ("ERROR: GeneTreeProbability.calculateGTDistribution(...) didn't return exactly one probability. Returning -1 to signal error.");
	    return (-1.0);
	}

        return (probList.get(0));
    }

    /**
     * Convert an EvoTree into a PhyloNet Tree object.
     */
    protected Tree convertPHMMTreeToPhyloNetTree (EvoTree etree) {
	NewickReader nr = new NewickReader(new StringReader(etree.toNewickString(true, true)));
	// kliu - change this over to a String data member
	// don't really use it anyways
	STITree<Double> newtr = new STITree<Double>(true);
	try {
	    nr.readTree(newtr);
	    return (newtr);
	}
	catch(Exception e) {
	    System.err.println(e);
	    e.printStackTrace();
	    return (null);
	}
    }

    /**
     * Convert an EvoTree into a PhyloNet Network object
     */
    protected Network<Double> convertPHMMTreeToPhyloNetNetwork (EvoTree etree) {
	ExNewickReader<Double> enr = new ExNewickReader<Double>(new StringReader(etree.toNewickString(true, true)));
	// kliu - change this over to a String data member
	// don't really use it anyways
	try {
	    Network<Double> network = enr.readNetwork();
	    return (network);
	}
	catch(IOException ioe) {
	    System.err.println(ioe);
	    ioe.printStackTrace();
	    return (null);
	}
    }
	
    /**
     * TODO. For now, just a warning.
     */
    public ObservationMap generate ()
    {	
	System.err.println ("ERROR: OpdfMap.generate() not implemented yet. Returning null.");
	return (null);
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
