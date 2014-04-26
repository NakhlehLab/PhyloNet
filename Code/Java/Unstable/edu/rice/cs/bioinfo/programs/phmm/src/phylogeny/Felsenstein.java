/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package phylogeny;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Hashtable;
import java.io.StringReader;
import optimize.CalculationCache;
import substitutionModel.GTRSubstitutionModel;
import substitutionModel.SubstitutionModel;
import substitutionModel.NucleotideAlphabet;
import substitutionModel.Alphabet;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;


/**
 * Simple Felsenstein calculator class.
 *
 * WARNING: no support for ambiguous nucleotides yet.
 * TODO: add this feature in? Nah.
 */

public class Felsenstein {
    protected static final String genes = NucleotideAlphabet.alphabet;

    protected SubstitutionModel substitutionModel;
    protected CalculationCache calculationCache;

    public Felsenstein (SubstitutionModel inSubstitutionModel,
			CalculationCache inCalculationCache) {
	this.substitutionModel = inSubstitutionModel;
	this.calculationCache = inCalculationCache;
    }

    // kliu - CACHE OPPORTUNITY #2
    // if gene tree branch lengths and substitution rates don't change, then re-use cached value if it exists.

    /**
     * On-the-fly version of getLikelihood() implemented for PhyloNet Trees
     * Calculates emission probability
     * P[O_t | g(s_i), b_{g(s_i)}, \theta ]
     *
     * See writeup for details.
     */
    public double getLikelihoodtree (Tree atree, ObservationMap column) {
        double result = 0;
	double[] baseFrequencies = substitutionModel.getStationaryProbabilities();
        ArrayList<Double> rootLikelihoods = getLikelihood(atree.getRoot(), column);
        for (int i = 0; i < 4; i++) {
            result += baseFrequencies[i] * rootLikelihoods.get(i);
        }

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
    protected double getPij(char i, char j, double t, TNode n) {
	Alphabet alphabet = NucleotideAlphabet.getClassInstance();

	// if cache entry exists, use it
	if (calculationCache.cacheSubstitutionProbabilityMatrix.containsKey(n)) {
	    double[][] pij = calculationCache.cacheSubstitutionProbabilityMatrix.get(n);
	    return (pij[alphabet.getObservationSymbolIndex(i)][alphabet.getObservationSymbolIndex(j)]);
	}

	// otherwise calculate and cache it
	double[][] pij = substitutionModel.calculateProbabilitiesFromRates(t);
	calculationCache.cacheSubstitutionProbabilityMatrix.put(n, pij);
	return (pij[alphabet.getObservationSymbolIndex(i)][alphabet.getObservationSymbolIndex(j)]);
    }


    /**
     * @return The node's likelihood arraylist
     */
    protected ArrayList<Double> getLikelihood(TNode aNode, ObservationMap column) {
        ArrayList<Double> likelihood = new ArrayList<Double>();

        if (aNode.isLeaf()) {

            //The observation for this leaf
            char obs = column.get(aNode.getName());

            // sets likelihood of leaf based on observation
            int index = -1;
            for (int i = 0; i < genes.length(); i++) {
                if (obs == genes.charAt(i)) {
		    index = i;
                }
            }

            // set likelihood for the observed observation as 1.0
            for (int i = 0; i < 4; i++) {
                if (i == index) {
		    likelihood.add(i, 1.0);
                }
                else likelihood.add(i, 0.0);
            }

            //System.out.println(" I am leaf : " + aNode.getName() + " and my likelihood array is : " + likelihood);
            return likelihood;
        }


        else {

	    Iterator<? extends TNode> childrenIterator= aNode.getChildren().iterator();
	    HashMap<TNode,ArrayList<Double>> childrenLikelihood = new HashMap<TNode,ArrayList<Double>>();
	    ArrayList<TNode> children = new ArrayList<TNode>();

	    //get all the children's likelihood arrays
	    while(childrenIterator.hasNext()) {
		TNode next = childrenIterator.next();
		children.add(next);
		childrenLikelihood.put(next, getLikelihood(next, column));
	    }

	    //Calculate likelihood for this node using Felsenstein's algorithm
	    for (int i = 0; i < genes.length(); i++) {
                double[] tempvalues = new double[aNode.getChildCount()];

                for (int j = 0; j < genes.length(); j++) {
                    for (int k = 0; k < tempvalues.length; k++) {
                        double childLikelihood = childrenLikelihood.get(children.get(k)).get(j);
			// last argument is for caching behavior
                        tempvalues[k] += childLikelihood * getPij(genes.charAt(i), genes.charAt(j), children.get(k).getParentDistance(), children.get(k));
                    }
                }

                double product = 1;
                for (int k = 0; k < tempvalues.length; k++) {
                    product *= tempvalues[k];
                }

                likelihood.add(i, product);
            }

            //System.out.println(" I am an internal node " + " and my likelihood array is : " + likelihood);

            return likelihood;
        }
    }

    // looks good
    protected static void test () {
	NewickReader nr = new NewickReader(new StringReader("((SA:1.0,SB:1.0)SAB:0.5,(SC:1.0,SD:1.0)SCD:0.5)root"));
	//NewickReader nr = new NewickReader(new StringReader("(((SA:1.0,SB:1.0)SAB:1.0,SC:1.0)SABC:0.5,SD:0.5)root"));
	// kliu - change this over to a String data member
	// don't really use it anyways
	STITree<Double> tree = new STITree<Double>(true);
	try {
	    nr.readTree(tree);
	}
	catch(Exception e) {
	    System.err.println(e);
	    e.printStackTrace();
	    return;
	}

	Map<String,Character> omap = new Hashtable<String,Character>();
	omap.put("SA", 'A');
	omap.put("SB", 'A');
	omap.put("SC", 'T');
	omap.put("SD", 'A');
	ObservationMap observationMap = new ObservationMap(omap);

	GTRSubstitutionModel gtrsm = new GTRSubstitutionModel();
	double[] rates = {1.0, 1.0, 1.0, 1.0, 1.0};
	double[] freqs = {0.25, 0.25, 0.25, 0.25};
	gtrsm.setSubstitutionRates(rates, freqs);

	Felsenstein fcalc = new Felsenstein(gtrsm, new CalculationCache());
	double result = fcalc.getLikelihoodtree(tree, observationMap);
	
	System.out.println ("Result: |" + result + "|");
    }

    public static void main (String[] args) {
	test();
    }

}


        // if (!i.equals(j)) {
        //     return (0.25 - 0.25 * Math.exp((-4.0 / 3.0) * u * t));
        // }
        // else
        //     return (0.25 + 0.75 * Math.exp((-4.0 / 3.0) * u * t));
