package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by yunyu on 4/22/15.
 */
public class Felsenstein {


    // kliu - CACHE OPPORTUNITY #2
    // if gene tree branch lengths and substitution rates don't change, then re-use cached value if it exists.

    /**
     * On-the-fly version of getLikelihood() implemented for PhyloNet Trees
     * Calculates emission probability
     * P[O_t | g(s_i), b_{g(s_i)}, \theta ]
     * <p>
     * See writeup for details.
     */
    public double getLikelihoodtree(Tree atree, Map<String, Character> map, SubstitutionModel model) {
        double result = 0;
        double[] baseFrequencies = {.1, .2, .3, .4};
        ArrayList<Double> rootLikelihoods = getLikelihood(atree.getRoot(), map, model);
        for (int i = 0; i < 4; i++) {
            result += baseFrequencies[i] * rootLikelihoods.get(i);
        }

        return result;
    }


    /**
     * Jukes-Cantor Model Calculation
     *
     * @param i - A Nucleotide Letter of type char {A, T, G, C}
     * @param j - A Nucleotide Letter of type char {A, T, G, C}
     * @param t - A Double - the time interval or branch length time
     *          return returns the Pij value, or the transition probability between two nucleotides given the branch length
     */
    protected double getPij(char i, char j, double t, SubstitutionModel model) {

        return model.getProbabilityMatrix(t).get(genes.indexOf(i), genes.indexOf(j));

    }


    public static final String genes = FelsensteinAlgorithm.nucleotides;

    /**
     * @return The node's likelihood ArrayList
     */
    protected ArrayList<Double> getLikelihood(TNode aNode, Map<String, Character> map, SubstitutionModel model) {
        ArrayList<Double> likelihood = new ArrayList<Double>();

        if (aNode.isLeaf()) {

            //The observation for this leaf
            char obs = map.get(aNode.getName());

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
                } else likelihood.add(i, 0.0);
            }

            //System.out.println(" I am leaf : " + aNode.getName() + " and my likelihood array is : " + likelihood);
            return likelihood;
        } else {

            Iterator<? extends TNode> childrenIterator = aNode.getChildren().iterator();
            HashMap<TNode, ArrayList<Double>> childrenLikelihood = new HashMap<TNode, ArrayList<Double>>();
            ArrayList<TNode> children = new ArrayList<TNode>();

            //get all the children's likelihood arrays
            while (childrenIterator.hasNext()) {
                TNode next = childrenIterator.next();
                children.add(next);
                childrenLikelihood.put(next, getLikelihood(next, map, model));
            }

            //Calculate likelihood for this node using Felsenstein's algorithm
            for (int i = 0; i < genes.length(); i++) {
                double[] tempvalues = new double[aNode.getChildCount()];

                for (int j = 0; j < genes.length(); j++) {
                    for (int k = 0; k < tempvalues.length; k++) {
                        double childLikelihood = childrenLikelihood.get(children.get(k)).get(j);
                        // last argument is for caching behavior
                        tempvalues[k] += childLikelihood * getPij( genes.charAt(j), genes.charAt(i), children.get(k).getParentDistance(), model);
                    }
                }

                double product = 1;
                for (double tempvalue : tempvalues) {
                    product *= tempvalue;
                }

                likelihood.add(i, product);
            }

            //System.out.println(" I am an internal node " + " and my likelihood array is : " + likelihood);

            return likelihood;
        }
    }
}
