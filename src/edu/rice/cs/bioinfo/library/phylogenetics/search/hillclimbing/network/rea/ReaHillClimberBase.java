package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAddition;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimberObservableBase;
import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Ref;

import java.util.Comparator;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/4/12
 * Time: 2:41 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReaHillClimberBase<G,N,E,S> extends HillClimberObservableBase<G, S>
{
    protected final ReticulateEdgeAddition<G,N,E> reaStrategy;

    private final boolean _performInitialNetworkValidation;

    private Map<N,Set<N>> _nodeToAncestors = null;

    protected Map<N,Set<N>> getNodeToAncestors()
    {
        return _nodeToAncestors;
    }

    private Func<Set<N>> _makeSet;

    public Func<Set<N>> getMakeSet()
    {
        return  _makeSet;
    }

    public void setMakeSet(Func<Set<N>> makeSet)
    {
        _makeSet = makeSet;
    }

    public ReaHillClimberBase(ReticulateEdgeAddition<G,N,E> reaStrategy, boolean performInitialNetworkValidation)
    {
        this.reaStrategy =  reaStrategy;
        _performInitialNetworkValidation = performInitialNetworkValidation;

        _makeSet = new Func<Set<N>>() {
            @Override
            public Set<N> execute() {
                return new HashSet<N>();
            }
        };
    }

    @Override
    protected boolean considerNeighborhood(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore) {
        final Ref<Boolean> sawBetterSolution = new Ref<Boolean>(false);
        final Ref<S> currentBestSeenSolutionScore = new Ref<S>(bestSeenSolutionScore);
        boolean shouldValidate = _performInitialNetworkValidation && this.getExaminationsCount() == 1;

       Func4<G,E,E,E,Boolean> rearrangementComputed = getRearrangementComputedListener(solution, getScore, scoreComparator, bestSeenSolutionScore,
                                                                                       getBetterSolution, newBestScore, sawBetterSolution, currentBestSeenSolutionScore);

       if(_nodeToAncestors == null)
       {
           _nodeToAncestors = reaStrategy.computeRearrangements(solution, shouldValidate, rearrangementComputed);
       }
       else
       {
           reaStrategy.computeRearrangements(solution, shouldValidate, rearrangementComputed, _nodeToAncestors);
       }
       return sawBetterSolution.get();
    }

    protected abstract Func4<G,E,E,E,Boolean> getRearrangementComputedListener(G solution, Func1<G, S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore,
                                                                               Ref<Func1<G, G>> getBetterSolution, Ref<S> newBestScore, Ref<Boolean> sawBetterSolution,
                                                                               Ref<S> currentBestSeenSolutionScore);
}
