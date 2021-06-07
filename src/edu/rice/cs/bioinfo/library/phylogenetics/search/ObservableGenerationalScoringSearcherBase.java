package edu.rice.cs.bioinfo.library.phylogenetics.search;

import edu.rice.cs.bioinfo.library.programming.Proc2;
import edu.rice.cs.bioinfo.library.programming.Proc4;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/8/12
 * Time: 11:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class ObservableGenerationalScoringSearcherBase<T,S> implements ObservableGenerationalScoringSearcher<T,S>
{
    private long _generationNumber = 0;

    protected long getGenerationNumber()
    {
        return _generationNumber;
    }

    protected void incrementGenerationNumber()
    {
        _generationNumber++;
    }

    protected void decrementGenerationNumber()
    {
        _generationNumber--;
    }

    private long _examinationsCount = 0;

    protected void incrementExaminations()
    {
        _examinationsCount++;
    }

    protected long getExaminationsCount()
    {
        return _examinationsCount;
    }

    private LinkedList<Proc4<T,S,Long, Long>> _betterSolutionFoundListeners = new LinkedList<Proc4<T,S,Long,Long>>();

    protected void fireBetterSolutionFoundEvent(T solution, S score)
    {
       for(Proc4<T,S,Long, Long> listener : _betterSolutionFoundListeners)
       {
           listener.execute(solution, score, getExaminationsCount(), getGenerationNumber());
       }
    }

    public void addBetterSolutionFoundListener(Proc4<T,S,Long,Long> listener)
    {
        if(!_betterSolutionFoundListeners.contains(listener))
        {
            _betterSolutionFoundListeners.add(listener);
        }
    }

    public boolean removeBetterSolutionFoundListener(Proc4<T,S,Long,Long> listener)
    {
        return _betterSolutionFoundListeners.remove(listener);
    }

    private LinkedList<Proc2<T,S>> _initialSolutionScoreComputedListeners = new LinkedList<Proc2<T, S>>();

    protected void fireInitialSolutionScoreComputedEvent(T solution, S score)
    {
        for(Proc2<T,S> listener : _initialSolutionScoreComputedListeners)
        {
            listener.execute(solution, score);
        }
    }

    public void addInitialSolutionScoreComputedListener(Proc2<T,S> listener)
    {
        if(!_initialSolutionScoreComputedListeners.contains(listener))
        {
            _initialSolutionScoreComputedListeners.add(listener);
        }
    }

    public boolean removeInitialSolutionScoreComputedListener(Proc2<T,S> listener)
    {
        return _initialSolutionScoreComputedListeners.remove(listener);
    }
}
