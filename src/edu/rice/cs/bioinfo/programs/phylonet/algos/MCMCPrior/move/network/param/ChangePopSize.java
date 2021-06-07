package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.param;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.List;

/**
 * Created by wendingqiao on 3/3/16.
 *
 * select a random node from all nodes in network
 * slide the popSize param in a fixed window
 * note: if new popSize is negative, take the absolute value
 */
public class ChangePopSize extends NetworkOperator {

    private double _windowSize = Utils._POP_SIZE_WINDOW_SIZE;

    private Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> _targetEdge;
    private Double _oldPopSize;

    public ChangePopSize(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        // reset
        _targetEdge = null;
        _oldPopSize = null;

        if(Utils._CONST_POP_SIZE) {
            _targetEdge = null;
            _oldPopSize = Double.valueOf(_network.getNetwork().getRoot().getRootPopSize());
        } else {
            List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = Networks.getAllEdges(_network.getNetwork());
            int rand = Randomizer.getRandomInt(edges.size() + 1);
            if(rand == edges.size()) {
                _targetEdge = null;
                _oldPopSize = Double.valueOf(_network.getNetwork().getRoot().getRootPopSize());
                if(_oldPopSize == Double.NaN) throw new RuntimeException("Invalid root pop size!!!!");
            } else {
                _targetEdge = edges.get(Randomizer.getRandomInt(edges.size()));
                _oldPopSize = Double.valueOf(_targetEdge.Item1.getParentSupport(_targetEdge.Item2));
            }
        }

        double newPopSize = Math.abs(_oldPopSize.doubleValue() + (Randomizer.getRandomDouble() - 0.5) * _windowSize);
        setPopSize(newPopSize);

        _violate = false;
        return 0;
    }

    @Override
    public void undo() {
        setPopSize(_oldPopSize.doubleValue());
    }

    @Override
    public String getName() {
        return "Change-PopSize";
    }

    @Override
    public void optimize(final double logAlpha) {
        _windowSize *= Math.exp(Utils.calcDelta(this, logAlpha));
    }

    private void setPopSize(double popSize) {
        if(_targetEdge == null) {
            _network.getNetwork().getRoot().setRootPopSize(popSize);
        } else {
            _targetEdge.Item1.setParentSupport(_targetEdge.Item2, popSize);
        }
    }

    // test
    public static void main(String[] args) {
        Utils._CONST_POP_SIZE = true;
        {
            UltrametricNetwork net = new UltrametricNetwork("((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0;");
            System.out.println(net.getNetwork().getRoot().getRootPopSize());
            int runs = 10000;
            int counter = 0;
            for (int i = 0; i < runs; i++) {
                NetworkOperator op = new ChangePopSize(net);
                double logHR = op.propose();
                if (logHR != Utils.INVALID_MOVE) counter++;
            }
            System.out.println(net.getNetwork().getRoot().getRootPopSize());
            System.out.println(net.getNetwork().toString());
            System.out.println(counter == runs);
            System.out.printf("%d out of %d\n", counter, runs);
        }
        {
            UltrametricNetwork net = new UltrametricNetwork(
                    "(((R:4.410551343749997:0.054539171485632476)I8#H2:0.23213428125000046:0.3157371615010589:0.5)I2#H1:1.0573143750000016:0.1729990140886931:0.5,((I8#H2:0.7336986562500014:0.029958091635864367:0.5,L:5.144249999999999:0.18851912138979207)I7:0.2707500000000005:0.36196379689505004,((((G:4.410551343749997:0.22427601255712748)I10#H3:0.23213428125000046:0.13013244851268:0.5,Q:4.642685624999998:0.3566908292108646)I4:0.24435187500000044:0.21470306296895045,I2#H1:0.24435187500000044:0.1686513888794018:0.5)I9:0.2572125000000005:0.24368977546344403,((I10#H3:0.23213428125000046:0.11856749121704276:0.5,C:4.642685624999998:0.0775115864818376)I6:0.24435187500000044:0.11298918485082364,A:4.887037499999998:0.47647587103261724)I5:0.2572125000000005:0.1686346046965678)I3:0.2707500000000005:0.3217493694068022)I1:0.28500000000000014:0.3356051886470172)I0;"
            );
            int runs = 10000;
            int test = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new ChangePopSize(net);
                double logHR = op.propose();
                op.undo();
                if(net.getNetwork().getRoot().getRootPopSize() != Utils._POP_SIZE_MEAN) {
                    test--;
                }
                test++;
            }
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
        }
    }
}
