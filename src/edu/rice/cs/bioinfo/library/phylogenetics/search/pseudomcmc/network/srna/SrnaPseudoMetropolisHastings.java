package edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc.network.srna;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAddition;
import edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc.PseudoMetropolisHastingsBase;
import edu.rice.cs.bioinfo.library.programming.DeepCopyable;
import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Ref;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/11/12
 * Time: 2:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class SrnaPseudoMetropolisHastings<G extends DeepCopyable<G>,N,E,S> extends PseudoMetropolisHastingsBase<G,S>
{
    protected final ReticulateEdgeAddition<G,N,E> reaStrategy;

    private final boolean _performInitialNetworkValidation;

     private Func<Set<N>> _makeSet;

    public Func<Set<N>> getMakeSet()
    {
        return  _makeSet;
    }

    public SrnaPseudoMetropolisHastings(ReticulateEdgeAddition<G,N,E> reaStrategy, boolean performInitialNetworkValidation)
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
    protected G considerNeighborhood(G currentSolution)
    {

        boolean shouldValidate = _performInitialNetworkValidation && this.getExaminationsCount() == 1;

        final Ref<E> acceptNetworkSourceEdge = new Ref<E>(null);
        final Ref<E> acceptNetworkDestinationEdge = new Ref<E>(null);
        final Ref<Boolean> nextNetworkAccepted = new Ref<Boolean>(false);


        Func4<G,E,E,E,Boolean> _rearrangementComputed = new Func4<G, E, E, E, Boolean>() {
           public Boolean execute(G neighbor, E sourceEdge, E destinationEdge, final E reticulateEdge)
           {
              if(SrnaPseudoMetropolisHastings.this.considerSolution(neighbor))
              {
                  return true;
              }
              else
              {
                  acceptNetworkSourceEdge.set(sourceEdge);
                  acceptNetworkDestinationEdge.set(destinationEdge);
                  nextNetworkAccepted.set(true);
                  return false;
              }

           }
       };

       Map<N,Set<N>> nodeToAncestors = reaStrategy.computeRearrangements(currentSolution, shouldValidate, _rearrangementComputed);

       while(!nextNetworkAccepted.get())
       {
            reaStrategy.computeRearrangements(currentSolution, false, _rearrangementComputed, nodeToAncestors);
       }

       reaStrategy.performRearrangement(currentSolution, false, acceptNetworkSourceEdge.get(), acceptNetworkDestinationEdge.get(), nodeToAncestors, _makeSet);

       return currentSolution;


    }
}
