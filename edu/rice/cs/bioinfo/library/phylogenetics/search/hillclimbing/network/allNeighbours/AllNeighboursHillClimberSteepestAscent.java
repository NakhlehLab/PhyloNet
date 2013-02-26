package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkWholeNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Ref;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/22/12
 * Time: 2:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class AllNeighboursHillClimberSteepestAscent <G extends Graph<N,E>,N,E,S> extends AllNeighboursHillClimberBase<G,N,E,S> {

    public AllNeighboursHillClimberSteepestAscent(NetworkWholeNeighbourhoodGenerator<G, N, E> generator) {
        super(generator);
    }

    @Override
    protected Func4<G,Integer,E,E,Boolean> getRearrangementComputedListener(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, S bestSeenSolutionScore,
                                                                      final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore,
                                                                      final Ref<Boolean> sawBetterSolution, final Ref<S> currentBestSeenSolutionScore)
    {
        return new Func4<G, Integer, E, E, Boolean>() {
            public Boolean execute(G neighbor, final Integer operation, final E edge1, final E edge2)
            {
                if (AllNeighboursHillClimberSteepestAscent.this.considerSolution(neighbor, getScore, scoreComparator, currentBestSeenSolutionScore.get(), newBestScore))
                {
                    sawBetterSolution.set(true);
                    currentBestSeenSolutionScore.set(newBestScore.get());
                    getBetterSolution.set(new Func1<G, G>() {
                        public G execute(G network) {
                            return _networkGenerator.performRearrangement(network,operation,edge1,edge2);
                        }
                    });
                }
                return getContinueSearch();
            }
        };
    }
}
