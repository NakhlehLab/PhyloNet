package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.edmondsAlgo;

import org.hamcrest.*;

import static ch.lambdaj.Lambda.filter;
import static org.hamcrest.Matchers.*;
import static org.hamcrest.Matchers.isIn;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/21/12
 * Time: 6:49 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class EdmondsAlgoGibbons85<V,E,W extends Comparable<W>>
{
    public class BranchingResult
    {
        public final Set<E> Branching;

        public final W BranchWeight;

        public BranchingResult(Set<E> branching, W branchWeight) {
            Branching = branching;
            BranchWeight = branchWeight;
        }
    }

    private class CycleRecord
    {
        public final Set<E> Cycle;

        public final E MinWeightEdgeInCycle;

        public final Set<E> ETildeEdgesOfCycle;

        public final V CycleVertex;

        CycleRecord(Set<E> cycle, E minWeightEdgeInCycle, Set<E> eTildeEdgesOfCycle, V cycleVertex) {


            if(!cycle.contains(minWeightEdgeInCycle))
            {
                throw new IllegalArgumentException();
            }

            if(!cycle.containsAll(eTildeEdgesOfCycle))
            {
                throw new IllegalArgumentException();
            }

            Cycle = cycle;
            MinWeightEdgeInCycle = minWeightEdgeInCycle;
            ETildeEdgesOfCycle = eTildeEdgesOfCycle;
            CycleVertex = cycleVertex;
        }
    }

    private final W _zero;

    public EdmondsAlgoGibbons85(W zero)
    {
        _zero = zero;
    }

    public BranchingResult findAMinimumBranching(Set<V> vertices, Set<E> directedEdges)
    {
        W sumOfAllPositiveWeightsAccum = _zero;
        for(E edge : directedEdges)
        {
            W weightOfEdge = getWeightOfEdge(edge);
            if(weightOfEdge.compareTo(_zero) > 0)
            {
                sumOfAllPositiveWeightsAccum = addWeight(sumOfAllPositiveWeightsAccum, weightOfEdge);
            }
        }

        Map<E,W> edgeToWeights = new HashMap<E, W>();
        for(E edge : directedEdges)
        {
            W weightOfEdge = getWeightOfEdge(edge);
            W adjustedWeight = subtractWeight(sumOfAllPositiveWeightsAccum, weightOfEdge);
            edgeToWeights.put(edge, adjustedWeight);
        }

        BranchingResult maxResult = findAMaximumBranchingHelp(vertices, directedEdges, edgeToWeights);

        W branchingWeightAccum = _zero;
        for(E edge : maxResult.Branching)
        {
            W edgeWeight = getWeightOfEdge(edge);
            branchingWeightAccum = addWeight(branchingWeightAccum, edgeWeight);

        }

        return new BranchingResult(maxResult.Branching, branchingWeightAccum);
    }

    public BranchingResult findAMaximumBranching(Set<V> vertices, Set<E> directedEdges)
    {
        Map<E,W> edgeToWeights = new HashMap<E, W>();
        for(E edge : directedEdges)
        {
            edgeToWeights.put(edge, getWeightOfEdge(edge));
        }

        return findAMaximumBranchingHelp(vertices, directedEdges, edgeToWeights);
    }

    private BranchingResult findAMaximumBranchingHelp(Set<V> vertices, Set<E> directedEdges, final Map<E,W> edgeToWeights)
    {
        Map<E,E> createdEdgeToPrevGenerationAnalogAccum = new HashMap<E, E>();

        HashSet<V> bv = new HashSet(); // the vertex bucket
        HashSet<E> be = new HashSet<E>();  // the edge bucket. Those edges provisionally chosen for branching.

        /*
         * Each discovered cycle requires the transformation of the problem into a new generation (that is a new graph).
         * Track these generational graphs.
         *
         * When the text mentions Ei that is edgeGenerations.get(i).
         */
        List<Set<V>> vertexGenerations = new LinkedList<Set<V>>(Arrays.asList(vertices)); // Vi in the text
        List<Set<E>> edgeGenerations = new LinkedList<Set<E>>(Arrays.asList(directedEdges)); // Ei in the text
        Map<Integer,CycleRecord> generationNumberToCycleAccum = new HashMap<Integer, CycleRecord>(); // Ci in the text.  First found cycle at 1.

        int i = 0; // the current generation number (starts at 0, with an increase for each discovered cycle).
        Set<V> vi = vertexGenerations.get(i); // the vertices for the ith generation.  Variable updated as i changes.
        Set<E> ei = edgeGenerations.get(i); // the edges for the ith generation.  Variable updated as i changes.

        while(!bv.containsAll(vi)) // line 3 in text
        {
            V v = filter(not(isIn(bv)), vi).get(0); // some vertex in the current generation but not in bv.  line 4 in the text
            bv.add(v);
            E e =  findAMaxEdgeWithDestination(v, ei, edgeToWeights); // line 6 in the text;
            if(e == null)
            {
                continue;
            }

            W weightOfe = edgeToWeights.get(e);

            if(weightOfe.compareTo(_zero) <= 0)  // line 7 in the text
                continue;

            HashSet<E> cycle = tryWouldIntroduceCycle(e, be);
            if(cycle != null) // line 8 in the text
            {
                i = i + 1; // cycle found, update generation number
                TransformResult transformResult = transformToNextGeneration(vi, ei, cycle, edgeGenerations.get(i-1), i, edgeToWeights, bv, be,
                                                                            createdEdgeToPrevGenerationAnalogAccum);
                vi = transformResult.Vi;
                vertexGenerations.add(vi);
                ei = transformResult.Ei;
                edgeGenerations.add(ei);
                generationNumberToCycleAccum.put(i, new CycleRecord(cycle, transformResult.MinWeightEdgeOfCycle, transformResult.ETildeEdgesOfCycle, transformResult.CycleVertex));
            }
            be.add(e); // line 12 in the text
        } // "goto 3" in the text

        // this is not in text.  My bookkeeping.
        Map<E,E> createdEdgeToPrevGenerationAnalog = createdEdgeToPrevGenerationAnalogAccum; // no longer accuming, rename variable
        Map<Integer,CycleRecord> generationNumberToCycle = generationNumberToCycleAccum;

        while(i!=0)
        {
            CycleRecord ithCycle = generationNumberToCycle.get(i);
            Set<E> ci = ithCycle.Cycle;
            V ui = ithCycle.CycleVertex;
            boolean isRootOfOutTree = isRootOfOutTreeIn(ui, be);

            for(E edge : new LinkedList<E>(be)) // avoid concurrent modification problems
            {
                V source = getSource(edge);
                V destination = getDestination(edge);
                if(ui.equals(source) || ui.equals(destination))
                {
                    be.remove(edge);
                    be.add(createdEdgeToPrevGenerationAnalog.get(edge));
                }
            }

            if(isRootOfOutTree)
            {
                for(E cycleEdge : ci )
                {
                    if(!cycleEdge.equals(ithCycle.MinWeightEdgeInCycle))
                        be.add(cycleEdge);
                }
            }
            else
            {
                for(E cycleEdge : ci )
                {
                    if(!ithCycle.ETildeEdgesOfCycle.contains(cycleEdge))
                        be.add(cycleEdge);
                }
            }
            i = i - 1;
        }

        W branchingWeightAccum = _zero;
        for(E edge : be)
        {
            W edgeWeight = edgeToWeights.get(edge);
            branchingWeightAccum = addWeight(branchingWeightAccum, edgeWeight);

        }

        return new BranchingResult(be, branchingWeightAccum);
    }

    private boolean isRootOfOutTreeIn(V considredRoot, HashSet<E> edges)
    {
        return isRootOfOutTreeInHelp(considredRoot, new HashSet<E>(edges));
    }

    private boolean isRootOfOutTreeInHelp(V considredRoot, HashSet<E> edges)
    {
        LinkedList<V> childAccum = new LinkedList<V>();
        Iterator<E> edgeElements = edges.iterator();
        while(edgeElements.hasNext())
        {
            E edge = edgeElements.next();
            V source = getSource(edge);
            V destination = getDestination(edge);

            if(considredRoot.equals(destination))
                return false;

            if(considredRoot.equals(source))
            {
                childAccum.add(destination);
                edgeElements.remove();
            }
        }

        for(V child : childAccum)
        {
            if(!isRootOfOutTreeIn(child, edges))
                return false;
        }

        return true;
    }



    class TransformResult { final Set<V> Vi; final Set<E> Ei; final E MinWeightEdgeOfCycle; final Set<E> ETildeEdgesOfCycle; final V CycleVertex;
        TransformResult(Set<V> vi, Set<E> ei, E minWeightEdgeOfCycle, Set<E> eTildeEdgesOfCycle, V cycleVertex) {
        Vi = vi;
        Ei = ei;
        MinWeightEdgeOfCycle = minWeightEdgeOfCycle;
        ETildeEdgesOfCycle = eTildeEdgesOfCycle;
        CycleVertex = cycleVertex;
        E MinWeightEdgeOfCycle = minWeightEdgeOfCycle;
    } }

    private TransformResult transformToNextGeneration(Set<V> vi, Set<E> ei, Iterable<E> cycle, Set<E> edgesLastGeneration,
                                                      int i, Map<E, W> edgeToWeights, HashSet<V> bv, HashSet<E> be, Map<E, E> createdEdgeToPrevGenerationAnalogAccum)
    {
       /*
        * Determine what the new vi should be.
        */
        HashSet<V> viToBe = new HashSet<V>(vi);
        HashSet<V> cycleNodes = new HashSet<V>();
        for(E cycleEdge : cycle)   // ensure no cycle nodes in viToBe
        {
            V source = getSource(cycleEdge);
            V destination = getDestination(cycleEdge);
            viToBe.remove(source);
            viToBe.remove(destination);
            cycleNodes.add(source);
            cycleNodes.add(destination);
        }
        V ui = makeVertex();
        viToBe.add(ui);


        /*
        * Determine what the new ei should be.
        */
        HashSet<E> eiToBe = new HashSet<E>(ei);
        Iterator<E> eiElements = eiToBe.iterator();
        while(eiElements.hasNext()) // ensure no edges incident to any cycle node in eiToBe
        {
            E edge = eiElements.next();
            V source = getSource(edge);
            V destination = getDestination(edge);

            if(cycleNodes.contains(source) || cycleNodes.contains(destination))
            {
                eiElements.remove();
            }
        }

        Set<E> eTildeEdgesOfCycle = new HashSet<E>();
        E minWeightEdgeOfCycle = findMinWeightCycleEdge(cycle, edgeToWeights);
        W minCycleWeight = edgeToWeights.get(minWeightEdgeOfCycle);
        for(E lastGenEdge : edgesLastGeneration) // update eiToBe with transformed edges that used to be leading to or away from cycle nodes, but not both
        {
            V source = getSource(lastGenEdge);
            V destination = getDestination(lastGenEdge);
            E createdEdge; // if edge from last generation is to/from cycle replace with a new edge to/from ui. store new edge here.
            if(cycleNodes.contains(source) && !cycleNodes.contains(destination)) // update edge leaving cycle
            {
                E fromCycleVertex = makeEdge(ui, destination);
                W lastGenEdgeWeight = edgeToWeights.get(lastGenEdge);
                edgeToWeights.put(fromCycleVertex, lastGenEdgeWeight);
                eiToBe.add(fromCycleVertex);

                createdEdge = fromCycleVertex;

            }
            else if(!cycleNodes.contains(source) && cycleNodes.contains(destination))  // update edge arriving in cycle
            {
                E toCycleVertex = makeEdge(source, ui);
                eiToBe.add(toCycleVertex);

                W lastGenEdgeWeight = edgeToWeights.get(lastGenEdge);
                E eTilde = findCycleEdgeWithDestination(destination, cycle);
                eTildeEdgesOfCycle.add(eTilde);
                W eTildeWeight = getWeightOfEdge(eTilde);

                W toCycleVertexWeight = addWeight(subtractWeight(lastGenEdgeWeight, eTildeWeight), minCycleWeight);  // formula on page 44
                edgeToWeights.put(toCycleVertex, toCycleVertexWeight);

                createdEdge = toCycleVertex;

            }
            else
            {
                createdEdge = null;
            }

            if(createdEdge != null) // this is not reflected in the text.  needed bookkeeping for smoothness.
            {
                createdEdgeToPrevGenerationAnalogAccum.put(createdEdge, lastGenEdge);
            }

            if(be.contains(lastGenEdge) && createdEdge != null) // if an edge be holds was updated, update in be
            {
                be.remove(lastGenEdge);
                be.add(createdEdge);
            }
        }


        Iterator<V> bvElements = bv.iterator();
        while(bvElements.hasNext())  // no cycle nodes in bv
        {
            V bvElement = bvElements.next();
            if(cycleNodes.contains(bvElement))
                bvElements.remove();
        }

        for(E cycleEdge : cycle) // no cycle edges in be
        {
            be.remove(cycleEdge);
        }

        return new TransformResult(viToBe, eiToBe, minWeightEdgeOfCycle, eTildeEdgesOfCycle, ui);
    }

    private E findMinWeightCycleEdge(Iterable<E> cycleEdges, Map<E, W> edgeToWeights)
    {
        E minSeen = cycleEdges.iterator().next();
        W minSeenWeight = edgeToWeights.get(minSeen);

        for(E cycleEdge : cycleEdges)
        {
            W cycleEdgeWeight = edgeToWeights.get(cycleEdge);

            if(cycleEdgeWeight.compareTo(cycleEdgeWeight) < 0)
            {
                minSeen = cycleEdge;
                minSeenWeight = cycleEdgeWeight;
            }
        }

        return minSeen;
    }

    protected E findCycleEdgeWithDestination(V destination, Iterable<E> cycleEdges)
    {
        for(E cycleEdge : cycleEdges)
        {
            V edgeDestination = getDestination(cycleEdge);
            if(destination.equals(edgeDestination))
                return cycleEdge;

        }

        throw new IllegalArgumentException("Destination not found as destination in given edges.");
    }

    private E findAMaxEdgeWithDestination(final V destination, Set<E> edges, Map<E, W> edgeToWeights) {

        Iterable<E> edgesWithDestination = filter(new CustomTypeSafeMatcher<E>("with destination")
        {
            @Override
            protected boolean matchesSafely(E edge) {
                V edgeDestination = getDestination(edge);
                boolean match = destination.equals(edgeDestination);
                return match;
            }
        }, edges);



        E maxEdgeSeen;

        try
        {
            maxEdgeSeen = edgesWithDestination.iterator().next();
        }
        catch(NoSuchElementException e)
        {
            return null;
        }

        W maxWeightSeen = edgeToWeights.get(maxEdgeSeen);
        for(E edge : edgesWithDestination)
        {
            W edgeWeight = edgeToWeights.get(edge);
            if(maxWeightSeen.compareTo(edgeWeight) < 0)
            {
                maxEdgeSeen = edge;
                maxWeightSeen = edgeWeight;
            }
        }
        return maxEdgeSeen;
    }

    protected HashSet<E> tryWouldIntroduceCycle(E edge, HashSet<E> edges)
    {
       Map<V,List<V>> adjMap = new HashMap<V, List<V>>();
       HashSet<E> allEdges = new HashSet<E>(edges);
       allEdges.add(edge);

       for(E e : allEdges)
       {
           V source = getSource(e);
           V dest = getDestination(e);

           if(!adjMap.containsKey(source))
           {
               adjMap.put(source,  new LinkedList<V>());
           }

           adjMap.get(source).add(dest);

           if(!adjMap.containsKey(dest))
               adjMap.put(dest, new LinkedList<V>());
       }

       V destination = getDestination(edge);
       V source = getSource(edge);
       LinkedList<V> cycle = tryWouldIntroduceCycleHelp(destination, source, adjMap, new LinkedList<V>());


       if(cycle != null)
       {
           cycle.add(destination);
           HashSet<E> cycleEdges = new HashSet<E>();

           ListIterator<V> cycleNodes = cycle.listIterator();
           V cycleSource = cycleNodes.next();
           while(cycleNodes.hasNext())
           {
                V cycleDest = cycleNodes.next();

                for(E e : allEdges)
                {
                    if(cycleSource.equals(getSource(e)) && cycleDest.equals(getDestination(e)))
                    {
                        cycleEdges.add(e);
                        break;
                    }
                }

               cycleSource = cycleDest;
           }
           return cycleEdges;
       }
       else
       {
           return null;
       }
    }

    private LinkedList<V> tryWouldIntroduceCycleHelp(V currentNode, V stopNode,  Map<V,List<V>> adjMap, LinkedList<V> pathAccum)
    {
        pathAccum.add(currentNode);
        if(currentNode.equals(stopNode))
            return pathAccum;

        for(V descendant : adjMap.get(currentNode))
        {
            LinkedList<V> cycle = tryWouldIntroduceCycleHelp(descendant, stopNode, adjMap, new LinkedList<V>(pathAccum));
            if(cycle != null)
                return cycle;
        }

        return null;
    }

    protected abstract W addWeight(W w1, W w2);

    protected abstract W subtractWeight(W w1, W w2);

    protected abstract E makeEdge(V source, V destination);

    protected abstract V makeVertex();

    protected abstract W getWeightOfEdge(E edge);

    protected abstract V getDestination(E edge);

    protected abstract V getSource(E edge);
}
