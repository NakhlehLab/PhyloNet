package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/17/12
 * Time: 4:00 PM
 * To change this template use File | Settings | File Templates.
 */
public class FindAllPredecessors<G extends GraphReadOnly<N,E>, N,E> implements Func1<G, Map<N,Set<N>>>
{

    @Override
    public Map<N, Set<N>> execute(G network) {

        if(!network.isRooted())
        {
            throw new IllegalArgumentException("Given graph must be rooted");
        }

        HashMap<N, Set<N>> nodeToAncestors = new HashMap<N, Set<N>>();

        N root = new FindRoot<N>().execute(network);
        nodeToAncestors.put(root, new HashSet<N>());

        Queue<Tuple<N,N>> toExpand = new LinkedList<Tuple<N,N>>();

        for(E rootOutEdge : network.getIncidentEdges(root))
        {
            toExpand.add(network.getNodesOfEdge(rootOutEdge));
        }

        while(toExpand.size() > 0)
        {
            Tuple<N,N> newEdge = toExpand.remove();
            N sourceNode = newEdge.Item1;
            N destinationNode = newEdge.Item2;

            Set<N> destinationAncestors;
            if(nodeToAncestors.containsKey(destinationNode))
            {
                destinationAncestors = nodeToAncestors.get(destinationNode);
            }
            else
            {
                destinationAncestors = new HashSet<N>();

            }
            destinationAncestors.addAll(nodeToAncestors.get(sourceNode));
            destinationAncestors.add(sourceNode);
            nodeToAncestors.put(destinationNode, destinationAncestors);


            for(E edge : network.getIncidentEdges(destinationNode))
            {
                Tuple<N,N> nodesOfEdge = network.getNodesOfEdge(edge);

                if(nodesOfEdge.Item1.equals(destinationNode))
                {
                    toExpand.add(nodesOfEdge);
                }
            }
        }
        return nodeToAncestors;
    }
}
