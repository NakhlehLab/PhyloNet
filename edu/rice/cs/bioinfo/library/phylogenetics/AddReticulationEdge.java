package edu.rice.cs.bioinfo.library.phylogenetics;

import com.sun.org.apache.bcel.internal.generic.RETURN;
import edu.rice.cs.bioinfo.library.programming.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/27/12
 * Time: 2:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class AddReticulationEdge<N,E,T extends Comparable<T>> implements Proc4<Graph<N,E>, E, E, Boolean>
{
    private Func1<Graph<N,E>,N> _makeNode;

    private Func3<Graph<N,E>,N,N,E> _makeEdge;

    private Func2<GraphReadOnly<N,E>,E,T> _getBranchLength;

    private Func<T> _makeZero;

    private Func2<T,T,T> _add;

    private Func2<T,T,T> _subtract;

    private Func1<T,T> _half;

    private Proc3<GraphReadOnly<N,E>, E, T> _setBranchLength;

    public AddReticulationEdge(Func1<Graph<N,E>,N> makeNode, Func3<Graph<N,E>,N,N,E> makeEdge, Func2<GraphReadOnly<N,E>,E,T> getBranchLength,
                               Proc3<GraphReadOnly<N,E>, E, T> setBranchLength,
                               Func<T> makeZero, Func2<T,T,T> add,  Func2<T,T,T> subtract, Func1<T,T> half)
    {
        _makeNode = makeNode;
        _makeEdge = makeEdge;
        _getBranchLength = getBranchLength;
        _setBranchLength = setBranchLength;
        _makeZero = makeZero;
        _add = add;
        _subtract = subtract;
        _half = half;
    }

    public void execute(Graph<N,E> graph, E edge1, E edge2, Boolean performValidation)
    {

        N edge1Source = graph.getNodesOfEdge(edge1).Item1;
        N edge2Source = graph.getNodesOfEdge(edge2).Item1;
        N edge1Dest = graph.getNodesOfEdge(edge1).Item2;
        N edge2Dest = graph.getNodesOfEdge(edge2).Item2;

        N root = new FindRoot<N,E>().execute(graph);

        T edge1BranchLength = _getBranchLength.execute(graph, edge1);
        T edge2BranchLength = _getBranchLength.execute(graph, edge2);

        FindPathLength<N,E,T> findPathLength = new FindPathLength<N, E, T>();
        T pathLengthToEdge1Source = findPathLength.execute(graph, root, edge1Source , _getBranchLength, _makeZero, _add);
        T pathLengthToEdge2Source = findPathLength.execute(graph, root, edge2Source, _getBranchLength, _makeZero, _add);

        T pathLengthToEdge1Dest = _add.execute(pathLengthToEdge1Source, edge1BranchLength);
        T pathLengthToEdge2Dest = _add.execute(pathLengthToEdge2Source, edge2BranchLength);

        boolean edge1StartsFirst = pathLengthToEdge1Source.compareTo(pathLengthToEdge2Source) < 0;

        boolean overlap;
        T reticulationTime;
        if(edge1StartsFirst)
        {
            overlap = pathLengthToEdge1Dest.compareTo(pathLengthToEdge2Source) > 0;
            reticulationTime = computeReticulationTime(pathLengthToEdge2Source, pathLengthToEdge2Dest, pathLengthToEdge1Dest);
        }
        else
        {
            overlap = pathLengthToEdge2Dest.compareTo(pathLengthToEdge1Source) > 0;
            reticulationTime = computeReticulationTime(pathLengthToEdge1Source, pathLengthToEdge1Dest, pathLengthToEdge2Dest);
        }

        if(!overlap)
        {
            throw new IllegalArgumentException("No overlap in time between edge1 and edge2.");
        }


        N node1 = _makeNode.execute(graph);
        graph.addNode(node1);
        NodeInjector.NodeInjectorUndoAction<Graph<N,E>,N,E> undo1 = NodeInjector.injectNodeIntoEdge(graph, edge1, node1, _makeEdge, false);
        E newEdge1To = graph.getEdge(edge1Source, node1);
        E newEdge1From = graph.getEdge(node1, edge1Dest);
        T toNewNode1BranchLength = _subtract.execute(reticulationTime, pathLengthToEdge1Source);
        T fromNewNode1BranchLength = _subtract.execute(pathLengthToEdge1Dest, reticulationTime);
        _setBranchLength.execute(graph, newEdge1To, toNewNode1BranchLength);
        _setBranchLength.execute(graph, newEdge1From, fromNewNode1BranchLength);

        N node2 = _makeNode.execute(graph);
        graph.addNode(node2);
        NodeInjector.NodeInjectorUndoAction<Graph<N,E>,N,E> undo2 = NodeInjector.injectNodeIntoEdge(graph, edge2, node2, _makeEdge, false);
        E newEdge2To = graph.getEdge(edge2Source, node2);
        E newEdge2From = graph.getEdge(node2, edge2Dest);
        T toNewNode2BranchLength = _subtract.execute(reticulationTime, pathLengthToEdge2Source);
        T fromNewNode2BranchLength = _subtract.execute(pathLengthToEdge2Dest, reticulationTime);
        _setBranchLength.execute(graph, newEdge2To, toNewNode2BranchLength);
        _setBranchLength.execute(graph, newEdge2From, fromNewNode2BranchLength);

        E reticulateEdge = _makeEdge.execute(graph, node1, node2);
        graph.addEdge(reticulateEdge);
        _setBranchLength.execute(graph, reticulateEdge, _makeZero.execute());



        if(performValidation)
        {
            try
            {
                GraphValidator.assertValidGraph(graph);
            }
            catch(IllegalArgumentException e)
            {
                graph.removeEdge(reticulateEdge);
                undo2.undoInjection();
                graph.removeNode(node2);
                undo1.undoInjection();
                graph.removeNode(node1);
                throw e;
            }
        }

    }

    private T computeReticulationTime(T start, T end1, T end2)
    {
        T minEnd = min(end1, end2);
        T range = _subtract.execute(minEnd, start);
        T delta = _half.execute(range);
        return _add.execute(start, delta);
    }

    private T min(T a, T b)
    {
        return a.compareTo(b) < 0 ? a : b;
    }
}
