package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing;

import edu.rice.cs.bioinfo.library.programming.Func1;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 1:37 PM
 * To change this template use File | Settings | File Templates.
 */
public interface HillClimber<T,S>
{
    public HillClimbResult<T,S> search(T startSolution, Func1<T,S> getScore, Comparator<S> scoreComparator);

    public HillClimbResult<T,S> search(T startSolution, Func1<T,S> getScore, Comparator<S> scoreComparator,
                                       long maxExaminations);
}
