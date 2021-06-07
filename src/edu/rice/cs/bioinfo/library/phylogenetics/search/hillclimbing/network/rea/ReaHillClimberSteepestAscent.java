package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAddition;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Ref;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 1:02 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReaHillClimberSteepestAscent<G,N,E,S> extends ReaHillClimberBase<G,N,E,S>
{
    public ReaHillClimberSteepestAscent(ReticulateEdgeAddition<G, N, E> reaStrategy, boolean performInitialNetworkValidation) {
        super(reaStrategy, performInitialNetworkValidation);
    }

    @Override
    protected Func4<G,E,E,E,Boolean> getRearrangementComputedListener(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, S bestSeenSolutionScore,
                                                                      final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore,
                                                                      final Ref<Boolean> sawBetterSolution, final Ref<S> currentBestSeenSolutionScore)
    {
        return new Func4<G, E, E, E, Boolean>() {
           public Boolean execute(G neighbor, final E sourceEdge, final E destinationEdge, final E reticulateEdge)
           {
               if (ReaHillClimberSteepestAscent.this.considerSolution(neighbor, getScore, scoreComparator, currentBestSeenSolutionScore.get(), newBestScore))
               {
                   sawBetterSolution.set(true);
                   currentBestSeenSolutionScore.set(newBestScore.get());
                   getBetterSolution.set(new Func1<G, G>() {
                       public G execute(G network) {
                           return reaStrategy.performRearrangement(network, false, sourceEdge, destinationEdge, getNodeToAncestors(), getMakeSet());
                       }
                   });
               }

               return getContinueSearch();
           }
       };
    }
}
