package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model;

import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util.SpeciesComparision;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util.SpeciesTreeComputations;

public class ProcessOutput
{

    public ProcessOutput()
    {
    }


    public ProcessedOutputInformation getProcessedOutput(RawOutputInformation rawOutputInformation)
    {
        System.out.println("Note: Using a hard coded threshold: Check ProcessOutput.java");
        return getProcessedOutput(rawOutputInformation, 0.75);
    }

    public ProcessedOutputInformation getProcessedOutput(RawOutputInformation rawOutputInformation, double threshold)
    {
        ProcessedOutputInformation result = new ProcessedOutputInformation();

        result.mostLikelySpeciesTrees = SpeciesTreeComputations.getMostLikelyTrees(rawOutputInformation.posteriorProbabilityOfSpeciesTrees);
        result.thresholdedMostLikelySpeciesTrees = SpeciesTreeComputations.getMostLikelyTreesWithThreshold(rawOutputInformation.posteriorProbabilityOfSpeciesTrees,threshold);
        result.threshold = threshold;

        if (rawOutputInformation.actualSpeciesTrees != null)
        {
            result.percentErrorBetweenActualAndViterbi = SpeciesComparision.percentError(rawOutputInformation.viterbiSpeciesTrees, rawOutputInformation.actualSpeciesTrees);
            result.percentErrorBetweenActualAndMostLikely = SpeciesComparision.percentError(result.mostLikelySpeciesTrees, rawOutputInformation.actualSpeciesTrees);
            result.percentErrorBetweenActualAndThreshold = SpeciesComparision.percentErrorAndUnknown(result.thresholdedMostLikelySpeciesTrees,rawOutputInformation.actualSpeciesTrees);
            result.errorMatrixForActualAndViterbi = SpeciesComparision.errorMatrix(rawOutputInformation.viterbiSpeciesTrees,rawOutputInformation.actualSpeciesTrees,rawOutputInformation.formattedAlleleMappings.size());
            result.errorMatrixForActualAndThreshold = SpeciesComparision.errorMatrix(result.thresholdedMostLikelySpeciesTrees,rawOutputInformation.actualSpeciesTrees,rawOutputInformation.formattedAlleleMappings.size());
        }

        return result;
    }

}
