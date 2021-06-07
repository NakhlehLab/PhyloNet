package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAddition;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.KSteepestAscentBase;
import edu.rice.cs.bioinfo.library.programming.DeepCopyable;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func4;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/18/12
 * Time: 1:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class SrnaKSteepestAscent<G extends DeepCopyable<G>, N, E, S> extends KSteepestAscentBase<G,S>
{
    private final ReticulateEdgeAddition<G,N,E> _reaStrategy;

    private final boolean _performInitialNetworkValidation;



    public SrnaKSteepestAscent(int numPaths, Func2<S,S,S> scoreDelta, ReticulateEdgeAddition<G,N,E> reaStrategy, boolean performInitialNetworkValidation)
    {
        super(numPaths, scoreDelta);
        _reaStrategy =  reaStrategy;
        _performInitialNetworkValidation = performInitialNetworkValidation;
    }

    @Override
    protected void considerNeighborhood(G solution, final S solutionScore, final long generation)
    {
      Func4<G,E,E,E,Boolean> rearrangementComputed = new Func4<G, E, E, E, Boolean>() {
        @Override
        public Boolean execute(G neighbor, E arg2, E arg3, E arg4) {
            return considerSolution(neighbor, solutionScore, generation + 1);
        }
      };

      _reaStrategy.computeRearrangements(solution, _performInitialNetworkValidation && generation == 0, rearrangementComputed);
    }
}
