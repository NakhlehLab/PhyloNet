package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing;

import edu.rice.cs.bioinfo.library.programming.*;
import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 1:37 PM
 * To change this template use File | Settings | File Templates.
 */
public interface HillClimber<T>
{
    public <S> HillClimbResult<T,S> search(T startSolution, Func1<T,S> getScore, Comparator<S> scoreComparator,
                                           int maxIterations);
}
