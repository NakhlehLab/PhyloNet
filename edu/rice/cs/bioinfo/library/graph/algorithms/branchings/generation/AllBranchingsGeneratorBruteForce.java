package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.generation;

import edu.rice.cs.bioinfo.library.graph.algorithms.cycleDetection.DirectedCycleDetectionDFS;
import edu.rice.cs.bioinfo.library.math.discrete.Configurations;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/9/13
 * Time: 4:25 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class AllBranchingsGeneratorBruteForce<E> implements Iterable<Set<E>>
{
    private class InEdgePermutationPart
    {
        public final Object Destination;

        private final List<E> _inEdges;

        private int _currentInEdgeIndex = 0;

        private InEdgePermutationPart(Object destination, List<E> inEdges) {
            Destination = destination;
            _inEdges = inEdges;
        }

        public boolean canAdvance()
        {
            return _currentInEdgeIndex != _inEdges.size() -1;
        }

        public void advance()
        {
            if(canAdvance())
            {
                _currentInEdgeIndex++;
            }
            else
            {
                throw new IllegalStateException("Can not advance");
            }
        }

        public void reset()
        {
            _currentInEdgeIndex = 0;
        }

        public E getCurrentEdge()
        {
            return _inEdges.get(_currentInEdgeIndex);
        }
    }

    private final Set<E> _graphEdges;

    private final Configurations<E> _inEdgeConfigurations;

    public AllBranchingsGeneratorBruteForce(Set<E> graphEdges)
    {
        _graphEdges = graphEdges;
        HashMap<Object, Set<E>> nodeToInEdgesAndNull = new HashMap<Object, Set<E>>();

        for(E edge : _graphEdges)
        {
            Object source = getSource(edge);
            Object destination = getDestination(edge);

            if(!nodeToInEdgesAndNull.containsKey(source))
                nodeToInEdgesAndNull.put(source, new HashSet<E>(Arrays.asList((E)null)));

            if(!nodeToInEdgesAndNull.containsKey(destination))
                nodeToInEdgesAndNull.put(destination, new HashSet<E>(Arrays.asList((E)null)));


            nodeToInEdgesAndNull.get(destination).add(edge);
        }

        LinkedList<List<E>> inEdgeChoicesForNode = new LinkedList<List<E>>();
        for(Object destination : nodeToInEdgesAndNull.keySet())
        {
            inEdgeChoicesForNode.add(new LinkedList<E>(nodeToInEdgesAndNull.get(destination)));

        }

        _inEdgeConfigurations = new Configurations<E>(inEdgeChoicesForNode);
    }

    public Iterator<Set<E>> iterator()
    {
        return new Iterator<Set<E>>() {

            private final Iterator<List<E>> _inEdgeConfigurations = AllBranchingsGeneratorBruteForce.this._inEdgeConfigurations.iterator();

            private List<E> _next = findNext(_inEdgeConfigurations);

            public boolean hasNext() {
                return _next != null;
            }

            public Set<E> next() {
               try
               {
                    return new HashSet<E>(_next);
               }
               finally
               {
                   _next = findNext(_inEdgeConfigurations);
               }
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    private List<E> findNext(Iterator<List<E>> inEdgeConfigurations)
    {
        while(inEdgeConfigurations.hasNext())
        {
            List<E> potentialNext = inEdgeConfigurations.next();
            potentialNext.removeAll(Arrays.asList((E)null));

            if(!containsCycle(potentialNext))
                return potentialNext;
        }

        return null;
    }

    private boolean containsCycle(List<E> edges) {

        return new DirectedCycleDetectionDFS<E>()
        {
            @Override
            protected Object getDestination(E edge) {
                return AllBranchingsGeneratorBruteForce.this.getDestination(edge);
            }

            @Override
            protected Object getSource(E edge) {
                return AllBranchingsGeneratorBruteForce.this.getSource(edge);
            }
        }.containsCycle( new HashSet<E>(edges));
    }

    protected abstract Object getSource(E edge);

    protected abstract Object getDestination(E edge);


}
