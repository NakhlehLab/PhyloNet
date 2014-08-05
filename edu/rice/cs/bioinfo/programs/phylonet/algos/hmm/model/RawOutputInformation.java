package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model;

public class RawOutputInformation
{
    public int[] viterbiSpeciesTrees;
    public int[] viterbiGeneTrees;

    public double[][] posteriorProbabilityOfSpeciesTrees;

    public double[][] posteriorProbabilityOfGeneTrees;

    public double finalScore;
    public double[] probScoreData;

    public double[] bestParametersRaw;
    public HmmParameters.Data bestParameters;


    public String[] formattedSpeciesTrees;
    public String formattedHmmStates;
    public int[] formattedHmmStateCounts;
    public String formattedNetwork;

    public int[] actualSpeciesTrees;
    public double scoreOfActual;
}
