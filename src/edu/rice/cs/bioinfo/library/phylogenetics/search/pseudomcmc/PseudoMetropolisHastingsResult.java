package edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/11/12
 * Time: 2:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class PseudoMetropolisHastingsResult<T, S>
{
    public final T BestExaminedNetwork;

    public final S BestExaminedScore;

    public final long ExaminationsCount;

    public final long GenerationCount;

    public PseudoMetropolisHastingsResult(T bestExaminedNetwork, S bestExaminedScore, long examinationCount, long generationCount)
    {
        BestExaminedNetwork = bestExaminedNetwork;
        BestExaminedScore = bestExaminedScore;
        ExaminationsCount = examinationCount;
        GenerationCount = generationCount;
    }
}