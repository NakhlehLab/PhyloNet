package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings;

import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.BranchingResult;
import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.MaxBranchingSolver;
import edu.rice.cs.bioinfo.library.programming.Container;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
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
public abstract class EdmondsAlgoGibbons85<V,E,W extends Comparable<W>> implements MaxBranchingSolver<V,E,W>
{
    private class DirectedEdge implements Comparable<DirectedEdge>
    {
        public final V Source;

        public final V Destination;

        public final W Weight;

        private DirectedEdge(V source, V destination, W weight) {
            Source = source;
            Destination = destination;
            Weight = weight;
        }

        public int compareTo(DirectedEdge o) {
            return this.Weight.compareTo(o.Weight);
        }
    }

    private class CycleRecord
    {
        public final Set<DirectedEdge> Cycle;

        public final DirectedEdge MinWeightEdgeInCycle;

        public final V CycleVertex;

        CycleRecord(Set<DirectedEdge> cycle, DirectedEdge minWeightEdgeInCycle, V cycleVertex) {


            if(!cycle.contains(minWeightEdgeInCycle))
            {
                throw new IllegalArgumentException();
            }

            Cycle = cycle;
            MinWeightEdgeInCycle = minWeightEdgeInCycle;
            CycleVertex = cycleVertex;
        }
    }

    private final W _zero;


    public EdmondsAlgoGibbons85(W zero)
    {
        _zero = zero;
    }

    public BranchingResult<E,W> findAMaximumBranching(Set<V> vertices, Set<E> directedEdges)
    {
        HashMap<DirectedEdge,E> directedEdgeToInputEdge = new HashMap<DirectedEdge,E>();
        for(E edge : directedEdges)
        {
            directedEdgeToInputEdge.put(new DirectedEdge(getSource(edge), getDestination(edge), getWeightOfEdge(edge)), edge);
        }

        return findAMaximumBranchingHelp(vertices, directedEdgeToInputEdge.keySet(), directedEdgeToInputEdge);
    }

    private BranchingResult<E,W> findAMaximumBranchingHelp(Set<V> vertices, Set<DirectedEdge> directedEdges, Map<DirectedEdge,E> directedEdgeToInputEdge)
    {
        Map<DirectedEdge,DirectedEdge> createdEdgeToPrevGenerationAnalogAccum = new HashMap<DirectedEdge, DirectedEdge>();

        HashSet<V> bv = new HashSet(); // the vertex bucket
        HashSet<DirectedEdge> be = new HashSet<DirectedEdge>();  // the edge bucket. Those edges provisionally chosen for branching.

        /*
         * Each discovered cycle requires the transformation of the problem into a new generation (that is a new graph).
         * Track these generational graphs.
         *
         * When the text mentions Ei that is edgeGenerations.get(i).
         */
        List<Set<V>> vertexGenerations = new LinkedList<Set<V>>(); // Vi in the text
        vertexGenerations.add(vertices);
        List<List<DirectedEdge>> edgeGenerations = new LinkedList<List<DirectedEdge>>(); // Ei in the text
        edgeGenerations.add(new LinkedList<DirectedEdge>(directedEdges));
        Map<Integer,CycleRecord> generationNumberToCycleAccum = new HashMap<Integer, CycleRecord>(); // Ci in the text.  First found cycle at 1.

        int i = 0; // the current generation number (starts at 0, with an increase for each discovered cycle).
        Set<V> vi = vertexGenerations.get(i); // the vertices for the ith generation.  Variable updated as i changes.
        List<DirectedEdge> ei = edgeGenerations.get(i); // the edges for the ith generation.  Variable updated as i changes.

        while(!bv.containsAll(vi)) // line 3 in text
        {
            V v = filter(not(isIn(bv)), vi).get(0); // some vertex in the current generation but not in bv.  line 4 in the text
            bv.add(v);
            DirectedEdge e =  tryFindAMaxEdgeWithDestination(v, ei); // line 6 in the text;
            if(e == null)
            {
                continue;
            }

            if(e.Weight.compareTo(_zero) <= 0)  // line 7 in the text
                continue;

            HashSet<DirectedEdge> cycle = tryWouldIntroduceCycle(e, be);
            if(cycle != null) // line 8 in the text
            {
                i = i + 1; // cycle found, update generation number
                TransformResult transformResult = transformToNextGeneration(vi, ei, cycle, edgeGenerations.get(i-1), i, bv, be,
                                                                            createdEdgeToPrevGenerationAnalogAccum);
                vi = transformResult.Vi;
                vertexGenerations.add(vi);
                ei = transformResult.Ei;
                edgeGenerations.add(ei);
                generationNumberToCycleAccum.put(i, new CycleRecord(cycle, transformResult.MinWeightEdgeOfCycle, transformResult.CycleVertex));
            }
            else
            {
                be.add(e); // line 12 in the text
            }
        } // "goto 3" in the text

        // this is not in text.  My bookkeeping.
        Map<DirectedEdge,DirectedEdge> createdEdgeToPrevGenerationAnalog = createdEdgeToPrevGenerationAnalogAccum; // no longer accuming, rename variable
        Map<Integer,CycleRecord> generationNumberToCycle = generationNumberToCycleAccum;

        while(i!=0)
        {
            CycleRecord ithCycle = generationNumberToCycle.get(i);
            Set<DirectedEdge> ci = ithCycle.Cycle;
            V ui = ithCycle.CycleVertex;

            DirectedEdge edgeToUi = tryFindInEdge(ui, be);

            for(DirectedEdge edge : new LinkedList<DirectedEdge>(be)) // avoid concurrent modification problems
            {

                if(ui.equals(edge.Source) || ui.equals(edge.Destination))
                {
                    be.remove(edge);
                    be.add(createdEdgeToPrevGenerationAnalog.get(edge));
                }
            }

            if(edgeToUi == null)
            {
                for(DirectedEdge cycleEdge : ci )
                {
                    if(!cycleEdge.equals(ithCycle.MinWeightEdgeInCycle))
                        be.add(cycleEdge);
                }
            }
            else
            {
                DirectedEdge inEdgeToCycleNode = createdEdgeToPrevGenerationAnalog.get(edgeToUi);
                V eTildeDestNode = inEdgeToCycleNode.Destination;
                for(DirectedEdge cycleEdge : ci )
                {
                    if(!cycleEdge.Destination.equals(eTildeDestNode))
                        be.add(cycleEdge);
                }
            }
            i = i - 1;
        }

        W branchingWeightAccum = _zero;
        HashSet<E> beOutsideContainer = new HashSet<E>();
        for(DirectedEdge edge : be)
        {
            branchingWeightAccum = addWeight(branchingWeightAccum, edge.Weight);
            beOutsideContainer.add(directedEdgeToInputEdge.get(edge));

        }

        return new BranchingResult(beOutsideContainer, branchingWeightAccum);
    }

    private DirectedEdge tryFindInEdge(V node, HashSet<DirectedEdge> edges)
    {
        for(DirectedEdge edge : edges)
        {
            if(edge.Destination.equals(node))
            {
                return edge;
            }
        }

        return null;
    }



    class TransformResult { final Set<V> Vi; final List<DirectedEdge> Ei; final DirectedEdge MinWeightEdgeOfCycle; final V CycleVertex;
        TransformResult(Set<V> vi, List<DirectedEdge> ei, DirectedEdge minWeightEdgeOfCycle, V cycleVertex) {
        Vi = vi;
        Ei = ei;
        MinWeightEdgeOfCycle = minWeightEdgeOfCycle;
        CycleVertex = cycleVertex;
        DirectedEdge MinWeightEdgeOfCycle = minWeightEdgeOfCycle;
    } }

    private TransformResult transformToNextGeneration(Set<V> vi, List<DirectedEdge> ei, Iterable<DirectedEdge> cycle, List<DirectedEdge> edgesLastGeneration,
                                                      int i, HashSet<V> bv, HashSet<DirectedEdge> be, Map<DirectedEdge, DirectedEdge> createdEdgeToPrevGenerationAnalogAccum)
    {
       /*
        * Determine what the new vi should be.
        */
        HashSet<V> viToBe = new HashSet<V>(vi);
        HashSet<V> cycleNodes = new HashSet<V>();
        for(DirectedEdge cycleEdge : cycle)   // ensure no cycle nodes in viToBe
        {
            viToBe.remove(cycleEdge.Source);
            viToBe.remove(cycleEdge.Destination);
            cycleNodes.add(cycleEdge.Source);
            cycleNodes.add(cycleEdge.Destination);
        }
        V ui = makeVertex();
        viToBe.add(ui);


        /*
        * Determine what the new ei should be.
        */
        List<DirectedEdge> eiToBe = new LinkedList<DirectedEdge>(ei);
        Iterator<DirectedEdge> eiElements = eiToBe.iterator();
        while(eiElements.hasNext()) // ensure no edges incident to any cycle node in eiToBe
        {
            DirectedEdge edge = eiElements.next();

            if(cycleNodes.contains(edge.Source) || cycleNodes.contains(edge.Destination))
            {
                eiElements.remove();
            }
        }

        //Set<E> eTildeEdgesOfCycle = new HashSet<E>();
        DirectedEdge minWeightEdgeOfCycle = IterableHelp.mins(cycle).iterator().next();
        for(DirectedEdge lastGenEdge : edgesLastGeneration) // update eiToBe with transformed edges that used to be leading to or away from cycle nodes, but not both
        {
            DirectedEdge createdEdge; // if edge from last generation is to/from cycle replace with a new edge to/from ui. store new edge here.
            if(cycleNodes.contains(lastGenEdge.Source) && !cycleNodes.contains(lastGenEdge.Destination)) // update edge leaving cycle
            {
                DirectedEdge fromCycleVertex = new DirectedEdge(ui, lastGenEdge.Destination, lastGenEdge.Weight);
                eiToBe.add(fromCycleVertex);
                createdEdge = fromCycleVertex;

            }
            else if(!cycleNodes.contains(lastGenEdge.Source) && cycleNodes.contains(lastGenEdge.Destination))  // update edge arriving in cycle
            {

                DirectedEdge eTilde = findCycleEdgeWithDestination(lastGenEdge.Destination, cycle);
                W toCycleVertexWeight = addWeight(subtractWeight(lastGenEdge.Weight, eTilde.Weight), minWeightEdgeOfCycle.Weight);  // formula on page 44

                DirectedEdge toCycleVertex =  new DirectedEdge(lastGenEdge.Source, ui, toCycleVertexWeight);
                eiToBe.add(toCycleVertex);

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

        for(DirectedEdge cycleEdge : cycle) // no cycle edges in be
        {
            be.remove(cycleEdge);
        }

        return new TransformResult(viToBe, eiToBe, minWeightEdgeOfCycle, ui);
    }
       /*
    private DirectedEdge findMinWeightEdge(Iterable<DirectedEdge> edges)
    {
        DirectedEdge minSeen = edges.iterator().next();

        for(DirectedEdge edge : edges)
        {
            if(edge.Weight.compareTo(minSeen.Weight) < 0)
            {
                minSeen = edge;
            }
        }

        return minSeen;
    }    */



    protected DirectedEdge findCycleEdgeWithDestination(V destination, Iterable<DirectedEdge> cycleEdges)
    {
        for(DirectedEdge cycleEdge : cycleEdges)
        {
            if(destination.equals(cycleEdge.Destination))
                return cycleEdge;

        }

        throw new IllegalArgumentException("Destination not found as destination in given edges.");
    }

    private DirectedEdge tryFindAMaxEdgeWithDestination(final V destination, List<DirectedEdge> edges) {

        Iterable<DirectedEdge> edgesWithDestination = filter(new CustomTypeSafeMatcher<DirectedEdge>("with destination")
        {
            @Override
            protected boolean matchesSafely(DirectedEdge edge) {
                boolean match = destination.equals(edge.Destination);
                return match;
            }
        }, edges);



        DirectedEdge maxEdgeSeen;

        try
        {
            maxEdgeSeen = edgesWithDestination.iterator().next();
        }
        catch(NoSuchElementException e)
        {
            return null;
        }

        for(DirectedEdge edge : edgesWithDestination)
        {
            if(maxEdgeSeen.Weight.compareTo(edge.Weight) < 0)
            {
                maxEdgeSeen = edge;
            }
        }
        return maxEdgeSeen;
    }

    protected HashSet<DirectedEdge> tryWouldIntroduceCycle(DirectedEdge edge, HashSet<DirectedEdge> edges)
    {
       Map<V,List<V>> adjMap = new HashMap<V, List<V>>();
       HashSet<DirectedEdge> allEdges = new HashSet<DirectedEdge>(edges);
       allEdges.add(edge);

       for(DirectedEdge e : allEdges)
       {
           if(!adjMap.containsKey(e.Source))
           {
               adjMap.put(e.Source,  new LinkedList<V>());
           }

           adjMap.get(e.Source).add(e.Destination);

           if(!adjMap.containsKey(e.Destination))
               adjMap.put(e.Destination, new LinkedList<V>());
       }


       LinkedList<V> cycle = tryWouldIntroduceCycleHelp(edge.Destination, edge.Source, adjMap, new LinkedList<V>());


        if(cycle != null)
        {
            cycle.add(edge.Destination);
            HashSet<DirectedEdge> cycleEdges = new HashSet<DirectedEdge>();

            ListIterator<V> cycleNodes = cycle.listIterator();
            V cycleSource = cycleNodes.next();
            while(cycleNodes.hasNext())
            {
                V cycleDest = cycleNodes.next();

                boolean foundMatch = false;
                for(DirectedEdge e : allEdges)
                {
                    if(cycleSource.equals(e.Source) && cycleDest.equals(e.Destination))
                    {
                        cycleEdges.add(e);
                        foundMatch = true;
                        break;
                    }
                }
               if(!foundMatch)
               {
                   throw new RuntimeException("Could not find match");
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

 //   protected abstract E makeEdge(V source, V destination);

    protected abstract V makeVertex();

    protected abstract W getWeightOfEdge(E edge);

    protected abstract V getSource(E edge);

    protected abstract V getDestination(E edge);


}
