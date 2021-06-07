package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.gui.PlotLineGraph;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.HmmTreeUtils;

import java.util.Arrays;
import java.util.List;

public class SpeciesComparision
{

    public static double[][] errorMatrix(int[] inferredSpeciesTree, int[] actualSpeciesTrees, int numberOfSpecies)
    {
        double[][] result = new double[numberOfSpecies][numberOfSpecies+1];
        for (int loci = 0; loci < inferredSpeciesTree.length; loci++)
        {
            if (inferredSpeciesTree[loci] == -1)
                result[actualSpeciesTrees[loci]][numberOfSpecies]+=1;
            else
                result[actualSpeciesTrees[loci]][inferredSpeciesTree[loci]]+=1;
        }

        for (int typeOfSpecies=0; typeOfSpecies < numberOfSpecies; typeOfSpecies++)
        {


            double sum = 0;
            for (int inferredSpeciesType=0; inferredSpeciesType <= numberOfSpecies; inferredSpeciesType++)
            {
                sum+=result[typeOfSpecies][inferredSpeciesType];
            }


            for (int inferredSpeciesType=0; inferredSpeciesType <= numberOfSpecies; inferredSpeciesType++)
            {
                if (sum != 0)
                    result[typeOfSpecies][inferredSpeciesType]/=sum;
                else
                    result[typeOfSpecies][inferredSpeciesType] = -1;
            }
        }
        return result;
    }

    public static class ErrorAndUnknown
    {
        public final double error;
        public final double unknown;

        public ErrorAndUnknown(double error, double unknown)
        {
            this.error = error;
            this.unknown = unknown;
        }
    }

    public static ErrorAndUnknown percentErrorAndUnknown(int[] inferredSpecies, int[] actualSpecies)
    {

        if (inferredSpecies.length != actualSpecies.length)
            throw new RuntimeException("Can only compare arrays of the same size: " + inferredSpecies.length + " " + actualSpecies.length);

        int errorCount = 0;
        int knownCount = 0;
        for (int i = 0; i < inferredSpecies.length;i++)
        {
            if (inferredSpecies[i] == -1)
                continue;

            knownCount++;
            if (inferredSpecies[i] != actualSpecies[i])
            {
                errorCount++;
            }

        }

        double error =  knownCount == 0 ? 1.0 : (double)errorCount/knownCount;
        double known = (double)knownCount/inferredSpecies.length;
        double unknown = 1- known;

        if (Double.isInfinite(error) || Double.isNaN(error))
            throw new RuntimeException("Error is not finite " + error + " " + unknown);

        if (Double.isInfinite(unknown) || Double.isNaN(unknown))
            throw new RuntimeException("Unknown is not finite " + error + " " + unknown);

        return new ErrorAndUnknown(error,unknown);
    }

    public static double percentError(int[] inferredSpecies, int[] actualSpecies)
    {
        ErrorAndUnknown result = percentErrorAndUnknown(inferredSpecies,actualSpecies);

        if (result.unknown != 0)
            throw new RuntimeException("There are unknowns (-1)s in the inferred species trees.");

        return result.error;
    }

    public static PlotLineGraph getCumulativeErrorPlot(int [] inferredSpecies, int[] actualSpecies)
    {
        if (inferredSpecies.length != actualSpecies.length)
            throw new RuntimeException("Can only compare arrays of the same size: " + inferredSpecies.length + " " + actualSpecies.length);


        double percentDifferent = 0.0;
        double increment = 1.0/inferredSpecies.length;
        double[] lineGraphData = new double[inferredSpecies.length];


        for (int i = 0; i < inferredSpecies.length; i++) {
            if (inferredSpecies[i] != actualSpecies[i])
                percentDifferent += increment;

            lineGraphData[i] = percentDifferent;
        }

        PlotLineGraph CStatPlot = new PlotLineGraph("Cumalitive Error of Viterbi Spec Trees compared to Reality Spec Trees",
                "Loci of Sequence", "Fraction of Viterbi Sequence != Reality Sequence", lineGraphData);
        return CStatPlot;
    }

    public static int[] getActualSpeciesTrees(List<Tree> geneTrees, List<STITree<String>> sourceTrees)
    {
        int[] result = new int[geneTrees.size()];
        for (int i = 0; i < result.length; i++)
        {
            int actualTree = getSourceTree(geneTrees.get(i),sourceTrees);
            result[i] = actualTree;
        }

        return result;
    }

    private static int getSourceTree(Tree geneTree, List<STITree<String>> sourceTrees)
    {
        int bestSpeciesIndex = -1;
        double bestSpeciesScore = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < sourceTrees.size(); i++)
        {
            double score = computeGeneTreeProbWithLenghts((STITree<String>) HmmTreeUtils.flattenTree(sourceTrees.get(i)), geneTree);
            if (score > bestSpeciesScore)
            {
                bestSpeciesScore = score;
                bestSpeciesIndex = i;
            }
        }

        if (bestSpeciesIndex == -1)
            throw new RuntimeException("Unable to find the correct source species tree");

        return bestSpeciesIndex;
    }

    private static double computeGeneTreeProbWithLenghts(STITree<String> speciesTree, Tree geneTree)
    {

        Network<Integer> test = HmmNetworkUtils.toNetwork(speciesTree);

        double[] result = new double[1];

        GeneTreeWithBranchLengthProbabilityYF g = new GeneTreeWithBranchLengthProbabilityYF(test, Arrays.asList(geneTree),null);
        g.calculateGTDistribution(result);
        return result[0];
    }

}
