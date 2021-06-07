package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 8:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class HillClimbResult<T, S>
{
    public final T BestExaminedNetwork;

    public final S BestExaminedScore;

    public final long ExaminationsCount;

    public final long GenerationCount;

    public HillClimbResult(T bestExaminedNetwork, S bestExaminedScore, long examinationCount, long generationCount)
    {
        BestExaminedNetwork = bestExaminedNetwork;
        BestExaminedScore = bestExaminedScore;
        ExaminationsCount = examinationCount;
        GenerationCount = generationCount;
    }
}
