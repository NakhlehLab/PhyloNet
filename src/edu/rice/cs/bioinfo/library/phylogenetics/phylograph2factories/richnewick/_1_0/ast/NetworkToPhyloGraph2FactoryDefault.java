package edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.programming.Func1;

import java.math.BigDecimal;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/12
 * Time: 1:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkToPhyloGraph2FactoryDefault<D> implements GraphFactory<RNNode, PhyloEdge2<RNNode,D>, Network>
{
    private final Func1<BigDecimal,D> _makeD;

    public NetworkToPhyloGraph2FactoryDefault(Func1<BigDecimal,D> makeD)
    {
        _makeD = makeD;
    }

    public PhyloGraph2<RNNode,D> make(Network input) {

        return input.execute(new NetworkAlgo<PhyloGraph2<RNNode,D>, Object, RuntimeException>() {
            public PhyloGraph2<RNNode,D> forNetworkEmpty(NetworkEmpty network, Object input) throws RuntimeException {
                return new DirectedPhyloGraphDefault2<RNNode,D>();
            }

            public PhyloGraph2<RNNode,D> forNetworkNonEmpty(NetworkNonEmpty network, Object input) throws RuntimeException {

                PhyloGraphBuilder<D> builder;
                if(network.RootageQualifier.execute(new IsRooted(), null))
                {
                    builder = new PhyloGraphBuilder<D>(new DirectedPhyloGraphDefault2<RNNode, D>(), _makeD);
                }
                else
                {
                    builder = new PhyloGraphBuilder<D>(new UndirectedPhyloGraphDefault2<RNNode, D>(), _makeD);
                }

                DAGFactory.makeDAG(network, builder);
                return builder.Graph;

            }
        }, null);
    }
}
