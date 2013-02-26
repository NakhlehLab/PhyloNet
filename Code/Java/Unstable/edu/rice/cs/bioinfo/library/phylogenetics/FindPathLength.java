package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func6;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/3/12
 * Time: 4:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class FindPathLength<N,E,T> implements Func6<GraphReadOnly<N,E>, N, N, Func2<GraphReadOnly<N,E>, E, T>, Func<T>, Func2<T,T,T>, T> {

    public T execute(GraphReadOnly<N,E> graph, N start, N end, Func2<GraphReadOnly<N,E>, E, T> getBranchLength, Func<T> makeZero, Func2<T,T,T> add)
    {
        if(!graph.isRooted())
        {
            throw new IllegalArgumentException("Given graph must be directed.");
        }

        T zero = makeZero.execute();
        T length = findPathLength(graph, start, end, getBranchLength, zero, add);

        if(length == null)
        {
            throw new IllegalArgumentException("End is not reachable from start.");
        }

        return length;
    }

    private T findPathLength(GraphReadOnly<N, E> graph, N current, N end, Func2<GraphReadOnly<N,E>, E, T> getBranchLength, T pathLengthAccum, Func2<T, T, T> add) {

        if(current.equals(end))
        {
            return pathLengthAccum;
        }

        for(N suc : new GetDirectSuccessors<N,E>().execute(graph, current))
        {
            E edge = graph.getEdge(current, suc);
            T branchLength =  getBranchLength.execute(graph, edge);
            T sucPathLength = add.execute(pathLengthAccum, branchLength);

            T length = findPathLength(graph, suc, end, getBranchLength, sucPathLength, add);

            if(length != null)
            {
                return length;
            }
        }

        return null;

    }


}
