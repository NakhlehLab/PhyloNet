package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing;

import edu.rice.cs.bioinfo.library.programming.*;

import java.util.Comparator;
import java.util.LinkedList;
import java.util.PriorityQueue;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/18/12
 * Time: 1:37 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class KSteepestAscentBase<T extends DeepCopyable<T>, S> implements HillClimberObservable<T,S>
{
    class SearchRecord
    {
        public final  T Solution;

        public final S Score;

        public final S ParentScoreDelta;

        public final long Generation;

        public SearchRecord(T solution, S score, S parentScoreDelta, Long generation)
        {
           Solution = solution;
           Score = score;
           ParentScoreDelta = parentScoreDelta;
           Generation = generation;
        }
    }

    private PriorityQueue<SearchRecord> _candidatePaths;

    private final int _numPaths;

    private Long _maxExaminations = null;

    private long _numExaminations = 0;

    private  long _largestGeneration = 0;

    private Func1<T, S> _getScore;

    private Comparator<S> _scoreComparator;

    private Func2<S,S,S> _scoreDelta;

    private SearchRecord _worstVelocityNetwork;

    private SearchRecord _bestScoreNetwork;

    public KSteepestAscentBase(int numPaths, Func2<S,S,S> scoreDelta)
    {
        _numPaths = numPaths;
        _scoreDelta = scoreDelta;
    }

    private LinkedList<Proc4<T, S, Long, Long>> _betterFoundListeners = new LinkedList<Proc4<T, S, Long, Long>>();

    public void addBetterSolutionFoundListener(Proc4<T, S, Long, Long> listener) {
        _betterFoundListeners.add(listener);
    }

    public boolean removeBetterSolutionFoundListener(Proc4<T, S, Long, Long> listener) {
        return _betterFoundListeners.remove(listener);
    }


    private LinkedList<Proc2<T, S>> _initialListeners = new LinkedList<Proc2<T, S>>();

    public void addInitialSolutionScoreComputedListener(Proc2<T, S> listener) {
        _initialListeners.add(listener);
    }

    public boolean removeInitialSolutionScoreComputedListener(Proc2<T, S> listener) {
        return _initialListeners.remove(listener);
    }

    public HillClimbResult<T, S> search(T startSolution, Func1<T, S> getScore, Comparator<S> scoreComparator, long maxExaminations) {
       _maxExaminations = maxExaminations;
       return search(startSolution, getScore, scoreComparator);

    }

    public HillClimbResult<T, S> search(T startSolution, Func1<T, S> getScore, final Comparator<S> scoreComparator)
    {
        _getScore = getScore;
        _scoreComparator = scoreComparator;

        Comparator<SearchRecord> reverseComparator = new Comparator<SearchRecord>() {
            @Override
            public int compare(SearchRecord o1, SearchRecord o2) {
                return scoreComparator.compare(o2.ParentScoreDelta, o1.ParentScoreDelta);
            }
        };
        _candidatePaths = new PriorityQueue<SearchRecord>(_numPaths, reverseComparator);

        S initialScore = getScore.execute(startSolution);
        _numExaminations++;
        _bestScoreNetwork = new SearchRecord(startSolution, initialScore, null, 0l);

        for(Proc2<T, S> listener : _initialListeners)
        {
            listener.execute(startSolution, initialScore);
        }

        considerNeighborhood(startSolution, initialScore, 0);

        while(_maxExaminations == null ? true : _numExaminations < _maxExaminations)
        {
            SearchRecord bestVelocity = _candidatePaths.poll();
            considerNeighborhood(bestVelocity.Solution, bestVelocity.Score, bestVelocity.Generation);
        }

        return new HillClimbResult<T,S>(_bestScoreNetwork.Solution, _bestScoreNetwork.Score, _numExaminations, _largestGeneration);
    }

    protected abstract void considerNeighborhood(T solution, S solutionScore, long generation);

    protected boolean considerSolution(final T solution, S searchParentScore, final long generation)
    {
        if(_numExaminations >= _maxExaminations)
        {
            return false;
        }

        final S score = _getScore.execute(solution);
        _numExaminations++;
        final S scoreDelta = _scoreDelta.execute(score, searchParentScore);

        if(generation > _largestGeneration)
        {
            _largestGeneration = generation;
        }

        Func<SearchRecord> makeSearchRecord = new Func<SearchRecord>() {
            @Override
            public SearchRecord execute() {
               return new SearchRecord(solution.DeepCopy(), score, scoreDelta, generation);
            }
        };

        if(_scoreComparator.compare(score, _bestScoreNetwork.Score) == 1)
        {
            _bestScoreNetwork = makeSearchRecord.execute();
            for(Proc4<T, S, Long, Long> listener : _betterFoundListeners)
            {
                listener.execute(solution, score, _numExaminations, generation);
            }
        }

        if(_worstVelocityNetwork == null)
        {
            SearchRecord candidatePath = makeSearchRecord.execute();
           _candidatePaths.add(candidatePath);
           _worstVelocityNetwork = candidatePath;
            return true;
        }


        boolean noBetterThanWorst = _scoreComparator.compare(scoreDelta, _worstVelocityNetwork.ParentScoreDelta) != 1;


        if(_candidatePaths.size() < _numPaths)
        {
           SearchRecord candidatePath = makeSearchRecord.execute();
           _candidatePaths.add(candidatePath);
           if(noBetterThanWorst)
           {
               _worstVelocityNetwork = candidatePath;
           }
           return true;
        }
        else if(noBetterThanWorst)
        {
            return true;
        }
        else
        {
             SearchRecord candidatePath = makeSearchRecord.execute();
            _candidatePaths.remove(_worstVelocityNetwork);
            _candidatePaths.add(candidatePath);
            _worstVelocityNetwork = _candidatePaths.peek();

            for(SearchRecord element : _candidatePaths)
            {
                if(_scoreComparator.compare(element.ParentScoreDelta, _worstVelocityNetwork.ParentScoreDelta) == -1)
                {
                    _worstVelocityNetwork = element;
                }
            }

            return true;
        }

    }
}
