package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.math.BigDecimal;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/4/12
 * Time: 6:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class IsUltrametric<E>
{
    public class Result
    {
        public final boolean IsUltrametricWithinThreshold;

        public final BigDecimal Height;

        Result(boolean isUltrametricWithinThreshold, BigDecimal height) {
            IsUltrametricWithinThreshold = isUltrametricWithinThreshold;
            Height = height;
        }
    }

    private final Func1<E, BigDecimal> _getBranchLength;

    public IsUltrametric(Func1<E, BigDecimal> getBranchLength)
    {
        _getBranchLength = getBranchLength;
    }

    public <N> Result execute(GraphReadOnly<N, E> graph, BigDecimal ultrametricThreshold)
    {
        Map<N, BigDecimal> pathLengthToLeafs = new HashMap<N, BigDecimal>();

               LinkedList<N> workingList = new LinkedList<N>();

               for(N leaf : new GetLeafs<N>().execute(graph))
               {
                   pathLengthToLeafs.put(leaf, BigDecimal.ZERO);
                   workingList.add(leaf);
               }

                GetDirectPredecessors<N,E> getDirectPredecessors = new GetDirectPredecessors<N, E>();

               N root = null;
               while(!workingList.isEmpty())
               {
                   N destNode = workingList.remove();
                   BigDecimal lengthFromDestNodeToLeaf = pathLengthToLeafs.get(destNode);
                   int numParentsForDestNode = 0;

                   for( N parent : getDirectPredecessors.execute(graph, destNode))
                   {
                       numParentsForDestNode++;
                       E incomingEdge = graph.getEdge(parent, destNode);
                       Tuple<N,N> nodesOfEdge = graph.getNodesOfEdge(incomingEdge);
                       N sourceNode = nodesOfEdge.Item1;
                       BigDecimal foundPathLength = pathLengthToLeafs.get(sourceNode);

                       BigDecimal incomingEdgeBranchLength = _getBranchLength.execute(incomingEdge);
                       BigDecimal expectedPathLength =  lengthFromDestNodeToLeaf.add(incomingEdgeBranchLength);

                       if(foundPathLength == null)
                       {
                           pathLengthToLeafs.put(sourceNode, expectedPathLength);
                           workingList.add(sourceNode);
                       }
                       else if(foundPathLength.subtract(expectedPathLength).abs().compareTo(ultrametricThreshold) == 1)
                       {
                           return new Result(false, null);
                       }
                   }

                   if(numParentsForDestNode == 0)
                   {
                      root = destNode;
                   }

               }


        return new Result(true, pathLengthToLeafs.get(root));
    }
}
