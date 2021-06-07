package edu.rice.cs.bioinfo.library.graph.algorithms.cycleDetection;


import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/9/13
 * Time: 5:05 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class DirectedCycleDetectionDFS<E>
{
    public boolean containsCycle(Set<E> edges)
    {
        return tryFindCycle(edges) != null;
    }

    public List<E> tryFindCycle(Set<E> edges)
    {
        Map<Object,List<Object>> adjMap = new HashMap<Object, List<Object>>();

        for(E e : edges)
        {
            Object source = getSource(e);
            Object dest = getDestination(e);

            if(!adjMap.containsKey(source))
                adjMap.put(source, new LinkedList<Object>());

            if(!adjMap.containsKey(dest))
                adjMap.put(dest, new LinkedList<Object>());

            adjMap.get(source).add(dest);
        }

        Set<Object> knownDAGRootsAccum = new HashSet<Object>();

        for(E edge : edges)
        {
            Object source = getSource(edge);
            LinkedList<Object> cycle = tryWouldIntroduceCycleHelp(source, adjMap, new LinkedList<Object>(), knownDAGRootsAccum);
            if(cycle != null)
            {
                List<E> cycleEdges = new LinkedList<E>();

                ListIterator<Object> cycleNodes = cycle.listIterator();
                Object cycleSource = cycleNodes.next();
                while(cycleNodes.hasNext())
                {
                    Object cycleDest = cycleNodes.next();

                    for(E e : edges)
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
        }

        return null;
    }

    private LinkedList<Object> tryWouldIntroduceCycleHelp(Object node, Map<Object,List<Object>> adjMap, LinkedList<Object> pathAccum, Set<Object> knownDAGRootsAccum)
    {
        try
        {
            if(pathAccum.contains(node))
            {
                return pathAccum;
            }
        }
        finally
        {
            pathAccum.add(node);
        }

        if(knownDAGRootsAccum.contains(node))
            return null;

        for(Object descendant : adjMap.get(node))
        {
            LinkedList<Object> cycle = tryWouldIntroduceCycleHelp(descendant, adjMap, new LinkedList<Object>(pathAccum), knownDAGRootsAccum);
            if(cycle != null)
                return cycle;
        }

        knownDAGRootsAccum.add(node);
        return null;
    }

    protected abstract Object getDestination(E edge);

    protected abstract Object getSource(E edge);


}
