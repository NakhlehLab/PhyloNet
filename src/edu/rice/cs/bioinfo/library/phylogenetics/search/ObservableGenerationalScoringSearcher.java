package edu.rice.cs.bioinfo.library.phylogenetics.search;

import edu.rice.cs.bioinfo.library.programming.Proc2;
import edu.rice.cs.bioinfo.library.programming.Proc4;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/8/12
 * Time: 11:31 AM
 * To change this template use File | Settings | File Templates.
 */
public interface ObservableGenerationalScoringSearcher<T,S>
{

    public void addBetterSolutionFoundListener(Proc4<T,S,Long,Long> listener);

    public boolean removeBetterSolutionFoundListener(Proc4<T,S,Long,Long> listener);

    public void addInitialSolutionScoreComputedListener(Proc2<T,S> listener);

    public boolean removeInitialSolutionScoreComputedListener(Proc2<T,S> listener);
}
