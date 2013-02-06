package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings;

import edu.rice.cs.bioinfo.library.graph.algorithms.cycleDetection.DirectedCycleDetectionDFS;
import edu.rice.cs.bioinfo.library.programming.Predicate1;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;


import java.util.*;

/**
 * A {@link MaxBranchingSolver} based on Edmonds' algorithm as described by Gibbons85.
 *
 * @see "Gibbons, Alan. "Algorithmic Graph Theory." New York: Press Syndicate of the University of Cambridge, 1985. Print."
 */
public abstract class EdmondsAlgoGibbons85<E,W extends Comparable<W>> implements MaxBranchingSolver<E,W>
{
    /**
     * Represents a weighted edge in a digraph.
     */
    private class DirectedEdge implements Comparable<DirectedEdge>
    {
        public final Object Source;

        public final Object Destination;

        public final W Weight;

        private DirectedEdge(Object source, Object destination, W weight) {
            Source = source;
            Destination = destination;
            Weight = weight;
        }

        /**
         * Compares the edge to the given edge wrt weight.
         */
        public int compareTo(DirectedEdge o) {
            return this.Weight.compareTo(o.Weight);
        }
    }

    /**
     * Represents a cycle discovered by the Edmonds algo.
     */
    private class CycleRecord
    {
        /**
         * The edges of the cycle.
         */
        public final Set<DirectedEdge> CycleEdges;

        /**
         * An edge of minimum weight in CycleEdges.
         */
        public final DirectedEdge MinWeightEdgeInCycle;

        /**
         * The single object node that represents this entire cycle in the next generation graph.
         * Referred to as ui in Gibbons85 pseudocode.
         */
        public final Object CycleVertex;

        /**
         * Creates a new cycle record.
         *
         * @param cycleEdges the edges of the discovered cycle
         * @param minWeightEdgeInCycle an edge of minimum weight in CycleEdges
         * @param cycleVertex the single object node that represents this entire cycle in the next generation graph
         */
        CycleRecord(Set<DirectedEdge> cycleEdges, DirectedEdge minWeightEdgeInCycle, Object cycleVertex) {

            if(!cycleEdges.contains(minWeightEdgeInCycle))
            {
                throw new IllegalArgumentException("minWeightEdgeInCycle is not contained in cycleEdges.");
            }

            CycleEdges = cycleEdges;
            MinWeightEdgeInCycle = minWeightEdgeInCycle;
            CycleVertex = cycleVertex;
        }
    }

    /**
     * A representation of zero in type W.
     */
    private final W _zero;

    /**
     * Creates a new EdmondsAlgo solver.
     * @param zero a representation of zero in type W.
     */
    public EdmondsAlgoGibbons85(W zero)
    {
        _zero = zero;
    }

    // contract from interface
    public BranchingResult<E,W> findAMaximumBranching(Set<E> directedEdges)
    {
        Set<Object> incidentVertices = new HashSet<Object>();

        // make a DirectedEdge for each input edge
        HashMap<DirectedEdge,E> directedEdgeToInputEdge = new HashMap<DirectedEdge,E>();
        for(E edge : directedEdges)
        {
            Object source = getSource(edge);
            Object destination = getDestination(edge);
            incidentVertices.add(source);
            incidentVertices.add(destination);
            DirectedEdge edgeRecord = new DirectedEdge(source, destination, getWeightOfEdge(edge));
            directedEdgeToInputEdge.put(edgeRecord, edge);
        }

        return findAMaximumBranchingHelp(incidentVertices, directedEdgeToInputEdge.keySet(), directedEdgeToInputEdge);
    }

    private BranchingResult<E,W> findAMaximumBranchingHelp(Set<Object> vertices, Set<DirectedEdge> directedEdges, Map<DirectedEdge,E> directedEdgeToInputEdge)
    {
        // NOTE: This algo is taken from "Gibbons, Alan. "Algorithmic Graph Theory." New York: Press Syndicate of the University of Cambridge, 1985. Print." p. 42-49.
        // HOWEVER.  There is a bug in the pseudocode of the book.  Line 12 should be in an else clause wrt. to the if on line 8.  We fix this bug in our implementation.

        // Every time Edmonds algo discovers a cycle in graph Gn it creates a new graph Gn+1 where that cycle is collapsed into a single node u_n+1.
        // As such if an edge (x,y) pointed into the cycle from outside the cycle in G a new edge (x,u_n+1) is created in Gn+1.
        // This maps stores the relationship between those edge conversions in Gn and Gn+1.
        // In our example, an entry would be made (x,u_n+1) -> (x,y) when creating Gn+1.
        Map<DirectedEdge,DirectedEdge> createdEdgeToPrevGenerationAnalogAccum = new HashMap<DirectedEdge, DirectedEdge>();

        final HashSet<Object> bv = new HashSet<Object>(); // the vertex bucket. BV in Gibbons85
        HashSet<DirectedEdge> be = new HashSet<DirectedEdge>();  // the edge bucket. Those edges provisionally chosen for branching. BE in Gibbons85.

        /*
         * Each discovered cycle requires the transformation of the problem into a new generation (that is a new graph).
         * Track these generational graphs Gi = (Vi, Ei)
         */
        List<Set<?>> vertexGenerations = new LinkedList<Set<?>>(); // Vi in the text
        vertexGenerations.add(vertices);
        List<List<DirectedEdge>> edgeGenerations = new LinkedList<List<DirectedEdge>>(); // Ei in the text
        edgeGenerations.add(new LinkedList<DirectedEdge>(directedEdges));
        Map<Integer,CycleRecord> generationNumberToCycleAccum = new HashMap<Integer, CycleRecord>(); // Ci in the text. Ci = generationNumberToCycleAccum.get(i). First found cycle at 1. No value at 0.

        int i = 0; // the current generation number (starts at 0, with an increase for each discovered cycle). G_0 is input graph to this method.
        Set<?> vi = vertexGenerations.get(i); // the vertices for the ith generation.  Variable updated as i changes. Vi in Gibbons85.
        List<DirectedEdge> ei = edgeGenerations.get(i); // the edges for the ith generation.  Variable updated as i changes. Ei in Gibbons85.

        while(!bv.containsAll(vi)) // line 3 in Gibbons85.
        {
            Object v = IterableHelp.filterUnknown(vi, new Predicate1<Object>() { // some vertex in the current generation but not in bv.  line 4 in Gibbons85.
                public boolean execute(Object node) {
                    return !bv.contains(node);
                }
            }).iterator().next();
            bv.add(v); // line 5 in Gibbons85.
            DirectedEdge e =  tryFindAMaxEdgeWithDestination(v, ei); // line 6 in the Gibbons85.
            if(e == null) // if no edge exists with destination v.
            {
                continue;
            }

            if(e.Weight.compareTo(_zero) <= 0)  // exclude non-positive weight. line 7 in the Gibbons85.
                continue;

            Set<DirectedEdge> cycle = tryWouldIntroduceCycle(e, be); // would the union of e and BE form a cycle? null means no.  value is the cycle.
            if(cycle != null) // if cycle. line 8 in Gibbons85.
            {
                i = i + 1; // cycle found, update generation number

                /*
                 * Transform graph into new graph with the found cycle collapsed into a single node.
                 * Compute new vi, ei, Ci, and update createdEdgeToPrevGenerationAnalogAccum.
                 * Lines 10-11 in Gibbons85.
                 */
                TransformResult transformResult = transformToNextGeneration(vi, ei, cycle, edgeGenerations.get(i-1), i, bv, be,
                                                                            createdEdgeToPrevGenerationAnalogAccum);
                vi = transformResult.Vi;  // update vi to reflect the vertices of the new graph generation
                vertexGenerations.add(vi);
                ei = transformResult.Ei; // update ei to reflect the edges of the new graph generation
                edgeGenerations.add(ei);
                CycleRecord ci = new CycleRecord(cycle, transformResult.MinWeightEdgeOfCycle, transformResult.CycleVertex);
                generationNumberToCycleAccum.put(i, ci);
            }
            else // NOTE: This else does not appear in Gibbons85.  The text is in error.  Algo doesn't work unless this is an else.
            {
                be.add(e); // line 12 in the text
            }
        } // "goto 3" in the text

        // this is not in text.  My bookkeeping.
        Map<DirectedEdge,DirectedEdge> createdEdgeToPrevGenerationAnalog = createdEdgeToPrevGenerationAnalogAccum; // no longer accuming, rename variable
        Map<Integer,CycleRecord> generationNumberToCycle = generationNumberToCycleAccum;

        while(i!=0) // line 14 in Gibbons85.
        {
            CycleRecord ithCycle = generationNumberToCycle.get(i);
            Set<DirectedEdge> ci = ithCycle.CycleEdges; // Ci in Gibbons85.
            Object ui = ithCycle.CycleVertex;

            DirectedEdge edgeToUi = tryFindInEdge(ui, be); // try to find any edge leading into ui in be.

            // "rename some edges in BE" line 15 in Gibbons85.
            // i.e. if a source or destination of an edge in be is ui, replace the edge with the corresponding edge from the previous generation.
            for(DirectedEdge edge : new LinkedList<DirectedEdge>(be)) // avoid concurrent modification problems removing elements from be directly
            {
                if(ui.equals(edge.Source) || ui.equals(edge.Destination))
                {
                    be.remove(edge);
                    be.add(createdEdgeToPrevGenerationAnalog.get(edge));
                }
            }

            if(edgeToUi == null) // "was ui a root of an out-tree in BE?" line 16 in Gibbons85.  Note BE always consists of trees.
            {                    // as such we only need to check that ui was a root, not that it is also a member of a tree vs dag or general digraph.
                // add all cycle edges into be except one of (arbitrary) minimum weight edges of the cycle.
                for(DirectedEdge cycleEdge : ci )
                {
                    if(cycleEdge != ithCycle.MinWeightEdgeInCycle) // line 17 in Gibbons85
                        be.add(cycleEdge);
                }
            }
            else
            {
                DirectedEdge inEdgeToCycle = createdEdgeToPrevGenerationAnalog.get(edgeToUi);
                Object eTildeDestNode = inEdgeToCycle.Destination; // the destination of the eTilde_i node in Gibbons85.

                // add all cycle edges to be except eTilde_i.
                for(DirectedEdge cycleEdge : ci )
                {
                    if(cycleEdge.Destination != eTildeDestNode) // is cycle edge eTilde_i?  line 18 in Gibbons85.
                        be.add(cycleEdge);
                }
            }
            i = i - 1;
        }

        // at this point the max branching has been found.  it is stored in be.

        W branchingWeightAccum = _zero; // the sum of the weights of the branching edges
        HashSet<E> branching = new HashSet<E>();
        // don't return branching typed as DirectedEdge, return as E.  Find the corresponding E for each DirectedEdge.
        for(DirectedEdge edge : be)
        {
            branchingWeightAccum = addWeight(branchingWeightAccum, edge.Weight);
            E branchingEdge = directedEdgeToInputEdge.get(edge);
            branching.add(branchingEdge);
        }

        return new BranchingResult<E,W>(branching, branchingWeightAccum);
    }

    private DirectedEdge tryFindInEdge(Object node, HashSet<DirectedEdge> edges)
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



    class TransformResult { final Set<Object> Vi; final List<DirectedEdge> Ei; final DirectedEdge MinWeightEdgeOfCycle; final Object CycleVertex;
        TransformResult(Set<Object> vi, List<DirectedEdge> ei, DirectedEdge minWeightEdgeOfCycle, Object cycleVertex) {
        Vi = vi;
        Ei = ei;
        MinWeightEdgeOfCycle = minWeightEdgeOfCycle;
        CycleVertex = cycleVertex;
        DirectedEdge MinWeightEdgeOfCycle = minWeightEdgeOfCycle;
    } }

    private TransformResult transformToNextGeneration(Set<?> vi, List<DirectedEdge> ei, Iterable<DirectedEdge> cycle, List<DirectedEdge> edgesLastGeneration,
                                                      int i, HashSet<Object> bv, HashSet<DirectedEdge> be, Map<DirectedEdge, DirectedEdge> createdEdgeToPrevGenerationAnalogAccum)
    {
       /*
        * Determine what the new vi should be.
        *
        * "As might be imagined Gi contains every vertex of Gi-1 except for those in Ci.  Vi also includes the new vertex ui." Gibbons 85.
        */
        HashSet<Object> viToBe = new HashSet<Object>(vi);
        HashSet<Object> cycleNodes = new HashSet<Object>();
        for(DirectedEdge cycleEdge : cycle)   // ensure no cycle nodes in viToBe
        {
            viToBe.remove(cycleEdge.Source);
            viToBe.remove(cycleEdge.Destination);
            cycleNodes.add(cycleEdge.Source);
            cycleNodes.add(cycleEdge.Destination);
        }
        Object ui = new Object();
        viToBe.add(ui);


       /*
        * Determine what the new ei should be.
        *
        * "Ei contains every edge of Ei-1 except those with one or more end-points in Ci.
        *  We also add to Ei new edges as follows.  For every edge (x,y) in Ei
        *  [that leads to or from the cycle but is not a cycle edge, update the endpoint with ui].
        */
        // remove cycle edges from eiToBe
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

        // update edges to/from cycle
        DirectedEdge minWeightEdgeOfCycle = IterableHelp.mins(cycle).iterator().next(); // arbitrary min weight edge of cycle
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

                DirectedEdge toCycle =  new DirectedEdge(lastGenEdge.Source, ui, toCycleVertexWeight);
                eiToBe.add(toCycle);

                createdEdge = toCycle;

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

        // line 11 in Gibbons85.
        Iterator<Object> bvElements = bv.iterator();
        while(bvElements.hasNext())  // no cycle nodes in bv
        {
            Object bvElement = bvElements.next();
            if(cycleNodes.contains(bvElement))
                bvElements.remove();
        }

        // line 11 in Gibbons85.
        for(DirectedEdge cycleEdge : cycle) // no cycle edges in be
        {
            be.remove(cycleEdge);
        }

        return new TransformResult(viToBe, eiToBe, minWeightEdgeOfCycle, ui);
    }

    private DirectedEdge findCycleEdgeWithDestination(Object destination, Iterable<DirectedEdge> cycleEdges)
    {
        for(DirectedEdge cycleEdge : cycleEdges)
        {
            if(destination == cycleEdge.Destination)
                return cycleEdge;

        }

        throw new IllegalArgumentException("Destination not found as destination in given edges.");
    }

    private DirectedEdge tryFindAMaxEdgeWithDestination(final Object destination, List<DirectedEdge> edges) {

        Collection<DirectedEdge> edgesWithDestination = IterableHelp.filter(edges, new Predicate1<DirectedEdge>() {
            public boolean execute(DirectedEdge edge) {
                boolean match = destination.equals(edge.Destination);
                return match;
            }
        });

        if(edgesWithDestination.size() > 0)
        {
            return IterableHelp.max(edgesWithDestination);
        }
        else
        {
            return null;
        }
    }

    private Set<DirectedEdge> tryWouldIntroduceCycle(DirectedEdge edge, HashSet<DirectedEdge> acyclicEdges)
    {
       Map<Object,List<Object>> adjMap = new HashMap<Object, List<Object>>();
       HashSet<DirectedEdge> allEdges = new HashSet<DirectedEdge>(acyclicEdges);
       allEdges.add(edge);

       List<DirectedEdge> cycle = new DirectedCycleDetectionDFS<DirectedEdge>()
       {
           @Override
           protected Object getDestination(DirectedEdge edge) {
               return edge.Destination;  //To change body of implemented methods use File | Settings | File Templates.
           }

           @Override
           protected Object getSource(DirectedEdge edge) {
               return edge.Source;  //To change body of implemented methods use File | Settings | File Templates.
           }
       }.tryFindCycle(allEdges);

       if(cycle != null)
       {
            return new HashSet<DirectedEdge>(cycle);
       }
       else
       {
           return null;
       }
    }


    protected abstract W addWeight(W w1, W w2);

    protected abstract W subtractWeight(W w1, W w2);

    protected abstract W getWeightOfEdge(E edge);

    protected abstract Object getSource(E edge);

    protected abstract Object getDestination(E edge);


}
