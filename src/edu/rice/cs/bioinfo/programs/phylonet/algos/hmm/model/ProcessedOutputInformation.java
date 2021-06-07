package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model;

import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util.SpeciesComparision;

public class ProcessedOutputInformation
{
    public double percentErrorBetweenActualAndViterbi;
    public double percentErrorBetweenActualAndMostLikely;

    public double[][] errorMatrixForActualAndViterbi;

    public SpeciesComparision.ErrorAndUnknown percentErrorBetweenActualAndThreshold;

    public double threshold;

    public int[] mostLikelySpeciesTrees;
    public int[] thresholdedMostLikelySpeciesTrees;

    public double[][] errorMatrixForActualAndThreshold;
}
