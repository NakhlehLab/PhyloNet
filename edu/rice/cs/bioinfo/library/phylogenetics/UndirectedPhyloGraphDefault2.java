package edu.rice.cs.bioinfo.library.phylogenetics;

import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.UndirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/25/12
 * Time: 2:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class UndirectedPhyloGraphDefault2<N,D> extends UndirectedGraphToGraphAdapter<N,PhyloEdge2<N,D>> implements PhyloGraph2<N, D>
{
    public UndirectedPhyloGraphDefault2()
    {
        super( new UndirectedSparseGraph<N, PhyloEdge2<N,D>>(), new Func1<PhyloEdge2<N, D>, Tuple<N, N>>() {
            public Tuple<N, N> execute(PhyloEdge2<N, D> edge) {
                return edge.NodesOfEdge;
            }
        });
    }
}
