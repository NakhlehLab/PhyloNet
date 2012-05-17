package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.*;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.*;
import edu.rice.cs.bioinfo.library.programming.*;

import java.util.Comparator;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 1:02 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReaHillClimber<G,N,E> extends HillClimberBase<G, Map<N,Set<N>>>
{
    private final ReticulateEdgeAddition<G,N,E> _reaStrategy;

    private final boolean _performInitialNetworkValidation;

    private Func<Set<N>> _makeSet;

    public void setMakeSet(Func<Set<N>> makeSet)
    {
        _makeSet = makeSet;
    }

    public ReaHillClimber(ReticulateEdgeAddition<G,N,E> reaStrategy, boolean performInitialNetworkValidation)
    {
        _reaStrategy =  reaStrategy;
        _performInitialNetworkValidation = performInitialNetworkValidation;

        setMakeSet(new Func<Set<N>>() {
            @Override
            public Set<N> execute() {
                return new HashSet<N>();
            }
        });
    }


    @Override
    protected <S> boolean considerNeighborhood(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,
                                               int iterationNumber, final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore, final Ref<Map<N,Set<N>>> nodeToAncestors) {
        final Ref<Boolean> sawBetterSolution = new Ref<Boolean>(false);
        final Ref<S> currentBestSeenSolutionScore = new Ref<S>(bestSeenSolutionScore);
        boolean shouldValidate = _performInitialNetworkValidation && iterationNumber == 1;

        Proc4<G,E,E,E> rearrangementComputed =  new Proc4<G, E, E, E>() {
           public void execute(G neighbor, final E sourceEdge, final E destinationEdge, final E reticulateEdge)
           {
               if (ReaHillClimber.this.considerSolution(neighbor, getScore, scoreComparator, currentBestSeenSolutionScore.get(), newBestScore))
               {
                   sawBetterSolution.set(true);
                   currentBestSeenSolutionScore.set(newBestScore.get());
                   getBetterSolution.set(new Func1<G, G>() {
                       public G execute(G network) {
                           return _reaStrategy.performRearrangement(network, false, sourceEdge, destinationEdge, nodeToAncestors.get(), _makeSet);
                       }
                   });
               }
           }
       };

       if(nodeToAncestors.get() == null)
       {
           Map<N,Set<N>> nodeToAncestorsValue = _reaStrategy.computeRearrangements(solution, shouldValidate, rearrangementComputed);
           nodeToAncestors.set(nodeToAncestorsValue);
       }
       else
       {
           _reaStrategy.computeRearrangements(solution, shouldValidate, rearrangementComputed, nodeToAncestors.get());
       }
       return sawBetterSolution.get();
    }
}
