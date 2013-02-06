package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.TransMapResult;

import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 2:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrerResult<E> implements TransMapResult<E>
{
    public final Set<Set<E>> Solutions;

    public final int GeneticDistanceTotal;

    public final int SilentColonizationTotal;

    public SnitkinTransMapInferrerResult(Set<Set<E>> solutions, int geneticDistanceTotal, int silentColonizationTotal)
    {
        GeneticDistanceTotal = geneticDistanceTotal;
        SilentColonizationTotal = silentColonizationTotal;
        Solutions = solutions;
    }

    public Set<Set<E>> getSolutions() {
        return Solutions;
    }
}
