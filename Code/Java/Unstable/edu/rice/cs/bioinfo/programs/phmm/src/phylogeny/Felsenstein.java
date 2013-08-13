package phylogeny;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

public class Felsenstein {


    static private String[] genes = {"A", "C", "G", "T"};


    // kliu - CACHE OPPORTUNITY #2
    // if gene tree branch lengths and substitution rates don't change, then re-use cached value if it exists.

     /**
     * On-the-fly version of getLikelihood() implemented for PhyloNet Trees
     * Calculates emission probability
     * P[O_t | g(s_i), b_{g(s_i)}, \theta ]
     *
     * See writeup for details.
     */
    static public double getLikelihoodtree (Tree atree, ObservationMap column, double baseSubRate) {

        double result = 0;
        ArrayList<Double> rootLikelihoods = getLikelihood(atree.getRoot(), column, baseSubRate);
        for (int i = 0; i < 4; i++) {
            result += 0.25 * rootLikelihoods.get(i);
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
    static private double getPij(String i, String j, double u, double t) {
        if (!i.equals(j)) {
            return (0.25 - 0.25 * Math.exp((-4.0 / 3.0) * u * t));
        }
        else
            return (0.25 + 0.75 * Math.exp((-4.0 / 3.0) * u * t));
    }


    /**
     * @return The node's likelihood arraylist
     */
    static private ArrayList<Double> getLikelihood(TNode aNode, ObservationMap column, double baseSubRate) {
        ArrayList<Double> likelihood = new ArrayList<Double>();

        if (aNode.isLeaf()) {

            //The observation for this leaf
            String obs = column.get(aNode.getName());

            // sets likelihood of leaf based on observation
            int index = -1;
            for (int i = 0; i < genes.length; i++) {
                if (obs.equals(genes[i])) {
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
                 childrenLikelihood.put(next, getLikelihood(next, column, baseSubRate));
             }

             //Calculate likelihood for this node using Felsenstein's algorithm
            for (int i = 0; i < genes.length; i++) {
                double[] tempvalues = new double[aNode.getChildCount()];

                for (int j = 0; j < genes.length; j++) {
                    for (int k = 0; k < tempvalues.length; k++) {
                        double childLikelihood = childrenLikelihood.get(children.get(k)).get(j);
                        tempvalues[k] += childLikelihood * getPij(genes[i], genes[j], baseSubRate, children.get(k).getParentDistance());
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


}
