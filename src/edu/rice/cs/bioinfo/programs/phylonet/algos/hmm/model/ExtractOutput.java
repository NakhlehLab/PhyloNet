package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import com.google.gson.Gson;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.CreateOptimumHMM;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util.SpeciesComparision;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util.StateProbabilityCalculator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class ExtractOutput
{

    HmmOptimizationFunction f;
    HmmParameters params;

    List<Tree> MSGeneTrees;
    List<STITree<String>> sourceTrees;

    public ExtractOutput(HmmOptimizationFunction the_f, List<Tree> MSGeneTrees, List<STITree<String>> sourceTrees)
    {
        f = the_f;
        params = f.getParams();
        this.MSGeneTrees = MSGeneTrees;
        this.sourceTrees = sourceTrees;
    }

    public RawOutputInformation getRawOutput(List<NucleotideObservation> testingObservations)
    {

        List<JahmmNucleotideObservation> JObservations = JahmmNucleotideObservation.wrapObservations(testingObservations);
        RawOutputInformation result = new RawOutputInformation();

        HmmParameters.Data hmmParameters = params.decode(f.getBestInput());
        Hmm<JahmmNucleotideObservation> hmm = f.getBestHmm();
        int[] viterbiStateSequence = hmm.mostLikelyStateSequence(JObservations);

        result.viterbiSpeciesTrees = toIntArray(HmmStateUtilities.getAlleleMappings(viterbiStateSequence, params));
        result.viterbiGeneTrees = toIntArray(HmmStateUtilities.getGeneTrees(viterbiStateSequence,params));

        result.posteriorProbabilityOfSpeciesTrees = getProbabilityOfSpeciesTreePerLoci(hmm, JObservations);

        result.posteriorProbabilityOfGeneTrees = getProbabilityOfGeneTreePerLoci(hmm);

        result.probScoreData = toDoubleArray(f.getProbScoreData());

        result.finalScore = f.getBestInputScore();
        result.bestParametersRaw = f.getBestInput();
        result.bestParameters = hmmParameters;

        result.formattedMulTree = f.getAlleleMappings(result.formattedAlleleMappings).toString();
        result.formattedHmmStates = getFormattedHmmStates(hmm);
        result.formattedHmmStateCounts =  getStateCounts(viterbiStateSequence,hmm.nbStates());

        if (MSGeneTrees != null)
            result.actualSpeciesTrees = getActualSpeciesTrees(MSGeneTrees);

        result.formattedNetwork = f.net.toString();

        if (result.actualSpeciesTrees != null)
        {
            double actualScore = probability(hmm, JObservations, result.actualSpeciesTrees);
            if (Double.isInfinite(actualScore))
                result.scoreOfActual = -1e200;
            else
                result.scoreOfActual = actualScore;
        }

        return result;
    }

    private int[] getStateCounts(int[] viterbiStateSequence, int numberOfStates)
    {
        int[] result = new int[numberOfStates];
        for (int state : viterbiStateSequence)
            result[state]++;

        return result;
    }

    private static double probability(Hmm<JahmmNucleotideObservation> hmm, List<JahmmNucleotideObservation> oseq, int[] sseq)
    {
        if (oseq.size() != sseq.length || oseq.isEmpty())
        {
            System.out.println("SSEQ Size = " +sseq.length);
            System.out.println("OSEQ Size = " +oseq.size());
            throw new IllegalArgumentException();
        }

        double probability = Math.log(hmm.getPi(sseq[0]));

        Iterator<JahmmNucleotideObservation> oseqIterator = oseq.iterator();

        for (int i = 0; i < sseq.length-1; i++)
            probability +=
                    Math.log(hmm.getOpdf(sseq[i]).probability(oseqIterator.next())) +
                            Math.log(hmm.getAij(sseq[i], sseq[i + 1]));

        return probability  + Math.log(hmm.getOpdf(sseq[sseq.length - 1]).
                probability(oseq.get(sseq.length - 1)));
    }


    private String[] getFormattedSpeciesTrees(List<STITree<String>> speciesTrees)
    {
        String[] result = new String[speciesTrees.size()];

        for (int i = 0 ; i < speciesTrees.size();i++)
            result[i] = speciesTrees.get(i).toNewickWD();

        return result;
    }

    private int[] getActualSpeciesTrees(List<Tree> actualMSGeneTrees)
    {
        return SpeciesComparision.getActualSpeciesTrees(actualMSGeneTrees, sourceTrees);

    }

    private String getFormattedHmmStates(Hmm<JahmmNucleotideObservation> hmm)
    {
        String result = "";
        for (int i = 0; i < hmm.nbStates();i++)
        {
            result += String.format("%d : %s\n", i, hmm.getOpdf(i));
        }
        return result;
    }

    private int[] toIntArray(List<Integer> intList) {
        int[] result = new int[intList.size()];

        for (int i = 0 ; i < intList.size();i++)
        {
            result[i] = intList.get(i);
        }

        return result;
    }

    private double[] toDoubleArray(List<Double> doubleList) {
        double[] result = new double[doubleList.size()];

        for (int i = 0 ; i < doubleList.size();i++)
        {
            result[i] = doubleList.get(i);
        }

        return result;
    }

    private double[][] getProbabilityOfSpeciesTreePerLoci(Hmm<JahmmNucleotideObservation> hmm, List<JahmmNucleotideObservation> JObservations) {
        StateProbabilityCalculator t = new StateProbabilityCalculator();

        double[][] gamma = t.getProbabilityOfStateAtTime(JObservations, hmm);
        double[][] prob = new double[params.getNumberOfAlleleMappings()][JObservations.size()];

        //Computation for a new plot.
        for (int i = 0; i < prob[0].length; i++) {
            for (int speciesTreeIndex = 0; speciesTreeIndex < params.getNumberOfAlleleMappings(); speciesTreeIndex++) {
                double[] val = new double[params.numberOfGeneTrees()];
                for (int geneTreeIndex = 0; geneTreeIndex < params.numberOfGeneTrees(); geneTreeIndex++) {
                    val[geneTreeIndex] = gamma[i][params.getStateIndex(geneTreeIndex, speciesTreeIndex)];
                }
                prob[speciesTreeIndex][i] = KahanSum(val);
            }
        }
        return prob;
    }

    private double[][] getProbabilityOfGeneTreePerLoci(Hmm<JahmmNucleotideObservation> hmm) {
        StateProbabilityCalculator t = new StateProbabilityCalculator();

        double[][] gamma = t.getProbabilityOfStateAtTime(f.getJahmmObservations(), hmm);
        double[][] prob = new double[params.numberOfGeneTrees()][f.getObservations().size()];

        //Computation for a new plot.
        for (int i = 0; i < prob[0].length; i++) {
            for (int geneTreeIndex = 0; geneTreeIndex < params.numberOfGeneTrees(); geneTreeIndex++) {
                double[] val = new double[params.getNumberOfAlleleMappings()];
                for (int speciesTreeIndex = 0; speciesTreeIndex < params.getNumberOfAlleleMappings(); speciesTreeIndex++) {
                    val[speciesTreeIndex] = gamma[i][params.getStateIndex(geneTreeIndex, speciesTreeIndex)];
                }
                prob[geneTreeIndex][i] = KahanSum(val);
            }
        }
        return prob;
    }

    private static double KahanSum(double[] input) {
        double sum = 0.0;
        double c = 0.0;
        for (double v : input) {
            double y = v - c;
            double t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        return sum;
    }


}
