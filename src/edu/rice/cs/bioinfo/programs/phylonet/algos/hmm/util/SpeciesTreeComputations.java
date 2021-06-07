package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.gui.DialogWithTabs;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.gui.PlotAreaGraph;


public class SpeciesTreeComputations {
    public static int[] getMostLikelyTrees(double[][] prob)
    {
        return getMostLikelyTreesWithThreshold(prob,0);
    }

    public static int[] getMostLikelyTreesWithThreshold(double[][] prob, double threshold)
    {
        int numTrees = prob.length;
        int numLoci = prob[0].length;
        int[] result = new int[numLoci];
        for (int locus = 0; locus < numLoci; locus++) {
            double maxProb = 0.0;
            int mostLikelyTree = -1;

            for (int tree = 0; tree < numTrees; tree++)
            {
                if(prob[tree][locus] > maxProb)
                {
                    maxProb = prob[tree][locus];
                    mostLikelyTree = tree;
                }
            }

            if (maxProb > threshold)
                result[locus] = mostLikelyTree;
            else
                result[locus] = -1;
        }

        return result;
    }

    public static double[][] getPostDecodeThresholdData(double[][] prob, int treeToThreshold, double threshold)
    {
        int numTrees = prob.length;
        int numLoci = prob[0].length;

        double[][] result = new double[numTrees][numLoci];

        //For each locus...
        for(int locus = 0; locus<numLoci; locus++)
        {
            if(prob[treeToThreshold][locus] >= threshold) //Set the result = the input @ this locus.
                for(int tree=0; tree<numTrees; tree++)
                    result[tree][locus] = prob[tree][locus];
            else //Set the result = 0 for the treeToThreshold, and = the input for others @ this locus.
            {
                for(int tree=0; tree<numTrees; tree++)
                {
                    if(tree == treeToThreshold)
                        result[tree][locus] = 0;
                    else
                        result[tree][locus] = prob[tree][locus];
                }
            }
        }
        return result;
    }

    public static double[][] getPostDecodeThresholdAllTreesData(double[][] prob, double threshold)
    {
        int numTrees = prob.length;
        int numLoci = prob[0].length;

        double[][] result = new double[numTrees][numLoci];

        //For each locus...
        for(int locus = 0; locus<numLoci; locus++)
        {
            //For each tree, erase data if below the threshold.
            for(int tree=0; tree<numTrees; tree++)
            {
                if(prob[tree][locus] >= threshold)
                    result[tree][locus] = prob[tree][locus];
                else
                    result[tree][locus] = 0;
            }
        }
        return result;
    }

    public static void addThresholdGraphPanel(DialogWithTabs dialogBox, double[][] prob, double threshold)
    {
        dialogBox.addPanel("Post-Decode w/Threshold", new PlotAreaGraph("Posterior Decoding Species Tree Probabilities"
        + "Threshold Analysis", getPostDecodeThresholdAllTreesData(prob, threshold)));
        /*for(int currentTree=0; currentTree<prob.length; currentTree++)
        {
            double[][] postDecodeThresholdData = getPostDecodeThresholdData(prob, currentTree, threshold);
            PlotAreaGraph thresholdGraph = new PlotAreaGraph("Posterior Decoding Species Tree Probabilities " +
                    "Threshold Analysis: Tree #" + currentTree, postDecodeThresholdData);
            dialogBox.addPanel("Post-Decode w/ Thresholds #" + currentTree, thresholdGraph);
        }*/
    }

}
