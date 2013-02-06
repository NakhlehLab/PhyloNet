package edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings;

import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.BranchingResult;
import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.EdmondsAlgoGibbons85;
import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.MaxBranchingSolver;
import edu.rice.cs.bioinfo.library.graph.algorithms.branchings.optimumBranchings.MaxBranchingSolverTest;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/7/13
 * Time: 7:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class EdmondsAlgoGibbons85Test extends MaxBranchingSolverTest
{


    @Override
    protected MaxBranchingSolver<String, Integer> makeSolverString(final Map<String, Integer> edgeToWeight) {
       return new EdmondsAlgoGibbons85<String,Integer>(0)
           {
               @Override
               protected Integer addWeight(Integer w1, Integer w2) {
                   return w1 + w2;
               }

               @Override
               protected Integer subtractWeight(Integer w1, Integer w2) {
                   return w1 - w2;
               }


               @Override
               protected Integer getWeightOfEdge(String edge) {
                   return edgeToWeight.get(edge);
               }

               @Override
               protected Character getDestination(String edge) {
                   return edge.charAt(1);  //To change body of implemented methods use File | Settings | File Templates.
               }

               @Override
               protected Character getSource(String edge) {
                   return edge.charAt(0);  //To change body of implemented methods use File | Settings | File Templates.
               }

       };

    }

    @Override
    protected MaxBranchingSolver<Tuple<Integer, Integer>, Integer> makeSolverTuple(final Map<Tuple<Integer, Integer>, Integer> edgeToWeight) {
        return new EdmondsAlgoGibbons85<Tuple<Integer, Integer>,Integer>(0)
        {
            @Override
            protected Integer addWeight(Integer w1, Integer w2) {
                return w1 + w2;
            }

            @Override
            protected Integer subtractWeight(Integer w1, Integer w2) {
                return w1 - w2;
            }

            @Override
            protected Integer getWeightOfEdge(Tuple<Integer, Integer> edge) {
                return edgeToWeight.get(edge);
            }

            @Override
            protected Integer getSource(Tuple<Integer, Integer> edge) {
                return edge.Item1;
            }

            @Override
            protected Integer getDestination(Tuple<Integer, Integer> edge) {
                return edge.Item2;
            }
        };
    }
}
