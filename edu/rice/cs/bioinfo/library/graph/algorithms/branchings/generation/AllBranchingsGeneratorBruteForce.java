package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.generation;

import edu.rice.cs.bioinfo.library.graph.algorithms.cycleDetection.DirectedCycleDetectionDFS;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/9/13
 * Time: 4:25 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class AllBranchingsGeneratorBruteForce<E>
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

    public AllBranchingsGeneratorBruteForce(Set<E> graphEdges)
    {
        _graphEdges = graphEdges;
    }


    public void generate()
    {
        Map<Object, Set<E>> nodeToInEdges = new HashMap<Object, Set<E>>();

        for(E edge : _graphEdges)
        {
            Object destination = getDestination(edge);

            if(!nodeToInEdges.containsKey(destination))
            {
                nodeToInEdges.put(destination, new HashSet<E>());
            }

            nodeToInEdges.get(destination).add(edge);
        }

        LinkedList<InEdgePermutationPart> permuntation = new LinkedList<InEdgePermutationPart>();

        for(Object node : nodeToInEdges.keySet())
        {
            Set<E> inEdges = nodeToInEdges.get(node);
            permuntation.add(new InEdgePermutationPart(node, new ArrayList<E>(inEdges)));
        }

        boolean moreToGenerate = true;

        while(moreToGenerate)
        {
            int numResets = 0;
            for(InEdgePermutationPart part : permuntation)
            {
                if(part.canAdvance())
                {
                    part.advance();
                    break;
                }
                else
                {
                    part.reset();
                    numResets++;
                }
            }

            moreToGenerate = numResets != permuntation.size();

            potentialBranchingGenerated(permuntation);
        }
    }

    private void potentialBranchingGenerated(LinkedList<InEdgePermutationPart> permutation)
    {
        Set<E> edges = new HashSet<E>();

        for(InEdgePermutationPart part : permutation)
        {
            edges.add(part.getCurrentEdge());
        }

        if(!containsCycle(edges))
        {
            branchingGenerated(edges);
        }
    }

    private boolean containsCycle(Set<E> edges) {

        final AllBranchingsGeneratorBruteForce outer = this;
        return new DirectedCycleDetectionDFS<E>()
        {
            @Override
            protected Object getDestination(E edge) {
                return outer.getDestination(edge);
            }

            @Override
            protected Object getSource(E edge) {
                return outer.getSource(edge);
            }
        }.containsCycle(edges);
    }


    protected abstract void branchingGenerated(Set<E> edges);

    protected abstract Object getDestination(E edge);

    protected abstract Object getSource(E edge);
}
