package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.Comparator;
import java.util.List;

/**
 * This class contains all of the output information for a given hill climbing attempt.
 */
public class FullOutputInformation implements Comparable<FullOutputInformation>
{
    public final RawOutputInformation rawOutputInformation;
    public final ProcessedOutputInformation processedOutputInformation;

    public FullOutputInformation(HmmOptimizationFunction the_f,List<Tree> MSTrees,List<STITree<String>> sourceTrees,List<NucleotideObservation> observations, double threshold)
    {
        rawOutputInformation = new ExtractOutput(the_f,MSTrees,sourceTrees).getRawOutput(observations);
        processedOutputInformation = new ProcessOutput().getProcessedOutput(rawOutputInformation, threshold);
    }

    /**
     * Creates a fulloutput information directly from a RawOutputInformation and ProcessedOutputInformation
     * @param rawOutputInformation The raw information for this run.
     * @param processedOutputInformation The processed information for this run.
     */
    public FullOutputInformation(RawOutputInformation rawOutputInformation, ProcessedOutputInformation processedOutputInformation)
    {
        this.rawOutputInformation = rawOutputInformation;
        this.processedOutputInformation = processedOutputInformation;
    }

    @Override
    public int compareTo(FullOutputInformation o)
    {
        return Double.compare(rawOutputInformation.finalScore,o.rawOutputInformation.finalScore);
    }
}