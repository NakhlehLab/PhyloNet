package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.*;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.*;
import edu.rice.cs.bioinfo.library.programming.*;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/13/12
 * Time: 1:02 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReaHillClimber<G,N,E> extends HillClimberBase<G>
{
    private final ReticulateEdgeAddition<G,N,E> _reaStrategy;

    public ReaHillClimber(ReticulateEdgeAddition<G,N,E> reaStrategy)
    {
        _reaStrategy =  reaStrategy;
    }


    @Override
    protected <S> boolean considerNeighborhood(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,
                                               final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore) {
        final Ref<Boolean> sawBetterSolution = new Ref<Boolean>(false);
       _reaStrategy.computeRearrangementsWithoutValidation(solution, new Proc4<G, E, E, E>() {
           public void execute(G neighbor, final E sourceEdge, final E destinationEdge, final E reticulateEdge) {

               if(ReaHillClimber.this.considerSolution(neighbor, getScore, scoreComparator, bestSeenSolutionScore, newBestScore))
               {
                   sawBetterSolution.set(true);
                   getBetterSolution.set(new Func1<G, G>()
                   {
                        public G execute(G network)
                        {
                            return _reaStrategy.perormRearrangementWithoutValidation(network, sourceEdge, destinationEdge);
                        }
                   });
               }
           }
       });
       return sawBetterSolution.get();
    }
}
