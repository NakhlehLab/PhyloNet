package edu.rice.cs.bioinfo.library.graph.algorithms.cycleDetection;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/16/13
 * Time: 1:17 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class CycleDetectorFromEdgesSets<E>
{
    public boolean containsCycle(Set<E> edges)
    {
        Hashtable<Object,Set<Object>> adjMap = new Hashtable<Object, Set<Object>>();

        for(E edge : edges)
        {
            Tuple<?,?> nodesOfEdge = getNodesOfEdge(edge);
            Object node1 = nodesOfEdge.Item1;
            Object node2 = nodesOfEdge.Item2;

            Set<Object> node1AdjSet = getSet(node1, adjMap);
            Set<Object> node2AdjSet = getSet(node2, adjMap);

            node1AdjSet.add(node2);
            node2AdjSet.add(node1);

        }

        if(adjMap.size() == 0)
            return false;

        HashSet<Object> potentialSearchRoots = new HashSet<Object>(adjMap.keySet());

        while(!potentialSearchRoots.isEmpty())
        {
            Object searchRoot = potentialSearchRoots.iterator().next();

            HashSet<Object> seenNodes = new HashSet<Object>();
            seenNodes.add(searchRoot);

            Stack<Object> toExplore = new Stack<Object>();
            toExplore.push(searchRoot);

            Stack<Object> searchParents = new Stack<Object>();
            searchParents.push(null);



            while(!toExplore.empty())
            {
                Object node = toExplore.pop();
                Object searchParent = searchParents.pop();
                potentialSearchRoots.remove(node);

                for(Object neighbor : adjMap.get(node))
                {
                    if(searchParent != null && searchParent == neighbor)
                        continue;

                    if(seenNodes.contains(neighbor))
                        return true;

                    seenNodes.add(neighbor);
                    toExplore.push(neighbor);
                    searchParents.push(node);
                }

            }
        }

        return false;
    }

    private Set<Object> getSet(Object node, Hashtable<Object,Set<Object>> adjMap)
    {
        Set<Object> set = adjMap.get(node);

        if(set == null)
        {
            set = new HashSet<Object>();
            adjMap.put(node, set);
        }

        return set;
    }

    protected abstract Tuple<?,?> getNodesOfEdge(E edge);
}
