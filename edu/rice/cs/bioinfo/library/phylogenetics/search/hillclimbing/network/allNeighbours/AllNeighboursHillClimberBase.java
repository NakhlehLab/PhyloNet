package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours.NetworkWholeNeighbourhoodGenerator;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimberObservableBase;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Ref;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/22/12
 * Time: 1:47 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class AllNeighboursHillClimberBase<G extends Graph<N,E>,N,E,S> extends HillClimberObservableBase<G, S> {
    protected int _diameterLimit;


    protected boolean considerNeighborhood(G solution, final Func1<G, S> getScore, final Comparator<S> scoreComparator, final S bestSeenSolutionScore,final Ref<Func1<G, G>> getBetterSolution, final Ref<S> newBestScore) {return true;}


    public HillClimbResult<G,S> search(G solution, Func1<G,S> getScore, Comparator<S> scoreComparator, Long maxExaminations, int maxGenerations, int diameterLimit)
    {
        _maxExaminations = maxExaminations;
        _maxGenerations = maxGenerations;
        _diameterLimit = diameterLimit;
        S score = getScore.execute(solution);
        this.fireInitialSolutionScoreComputedEvent(solution, score);
        incrementExaminations();
        return search(solution, getScore, scoreComparator, score);
    }


    public HillClimbResult<G,S> search(G solution, Func1<G,S> getScore, Comparator<S> scoreComparator, Long maxExaminations, int maxGenerations)
    {
        return search(solution, getScore, scoreComparator, maxExaminations, maxGenerations, 0);
    }



    protected abstract Func4<G,Integer,E,E,Boolean> getRearrangementComputedListener(G solution, Func1<G, S> getScore, Comparator<S> scoreComparator, S bestSeenSolutionScore,
                                                                               Ref<Func1<G, G>> getBetterSolution, Ref<S> newBestScore, Ref<Boolean> sawBetterSolution,
                                                                               Ref<S> currentBestSeenSolutionScore);

}
