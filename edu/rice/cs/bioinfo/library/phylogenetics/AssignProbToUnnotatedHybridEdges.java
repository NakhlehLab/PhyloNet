package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.math.BigDecimal;
import java.util.LinkedList;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/10/12
 * Time: 4:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class AssignProbToUnnotatedHybridEdges<N,E>  {

    private static Random _rand = new Random();

    public static void setRandom(Random newRand)
    {
        _rand = newRand;
    }

    public static Func2<Object,Object,BigDecimal> UNIFORM_RANDOM = new Func2<Object, Object, BigDecimal>() {
        public BigDecimal execute(Object graph, Object edge) {
           BigDecimal rand = new BigDecimal(_rand.nextDouble()); // 1 not possible. 0 is possible

           while(rand.compareTo(BigDecimal.ZERO) == 0) // hybrid edge with prob 0 (or 1) makes no sense.
           {
                rand = new BigDecimal(_rand.nextDouble());
           }

           return rand;
        }
    };

    public void execute(GraphReadOnly<N,E> graph, Func2<GraphReadOnly<N,E>,E,BigDecimal> getProb, Proc3<GraphReadOnly<N,E>, E, BigDecimal> setProb, Func2<GraphReadOnly<N,E>, E, Boolean> isEdgeProbUnset)
    {
        IsDestinationNode isDestNode = new IsDestinationNode();
        for(N node : graph.getNodes())
        {
            Iterable<E> edges = graph.getIncidentEdges(node);
            LinkedList<E> inEdges = new LinkedList<E>();

            for(E edge : edges)
            {
               if(isDestNode.execute(graph, node, edge))
               {
                  inEdges.add(edge);
               }
            }

            if(inEdges.size() == 2)
            {
                if(isEdgeProbUnset.execute(graph, inEdges.get(0)) && isEdgeProbUnset.execute(graph, inEdges.get(1)))
                {
                    BigDecimal prob = getProb.execute(graph, inEdges.get(0));

                    if(prob.compareTo(BigDecimal.ONE) == 1 || prob.compareTo(BigDecimal.ZERO) == 0)
                    {
                        throw new IllegalArgumentException("getProb generated a double outside the range [0,1].");
                    }

                    setProb.execute(graph, inEdges.get(0), prob);
                    setProb.execute(graph, inEdges.get(1), BigDecimal.ONE.subtract(prob));
                }
            }
        }
    }
}
