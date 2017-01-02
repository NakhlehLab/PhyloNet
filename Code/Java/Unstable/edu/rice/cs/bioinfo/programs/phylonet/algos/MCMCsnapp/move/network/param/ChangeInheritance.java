package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.param;

import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SpeciesNetPriorDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

import java.util.List;

/**
 * Created by wendingqiao on 3/5/16.
 *
 * select a random network node
 * propose a random rate within the sliding window
 * node: the rate should be in [0, 1], if rate > 1, rate = 2-rate; if rate < 0, rate = -rate
 */
public class ChangeInheritance extends NetworkOperator {

    private double _windowSize = 0.1;
    private Tuple3<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> _target;
    private Double _oldRate;

    public ChangeInheritance(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        // reset
        _target = null;
        _oldRate = null;

        List<NetNode<NetNodeInfo>> networkNodes = IterableHelp.toList(_network.getNetwork().getNetworkNodes());
        NetNode<NetNodeInfo> target = networkNodes.get(Randomizer.getRandomInt(networkNodes.size()));
        List<NetNode<NetNodeInfo>> parents = IterableHelp.toList(target.getParents());
        if(parents.size() != 2) {
            throw new RuntimeException("parents size() = " + parents.size() + "\n" + _network.getNetwork().toString());
        }
        if( Math.abs(target.getParentProbability(parents.get(0)) + target.getParentProbability(parents.get(1)) - 1.0) > 0.000001) {
            throw new RuntimeException("parents probability != 1.0 " + target.getName() + "\n" + _network.getNetwork().toString());
        }
        _target = new Tuple3<>(target, parents.get(0), parents.get(1));
        _oldRate = Double.valueOf(target.getParentProbability(_target.Item2));
        double newRate = _oldRate.doubleValue() + (Randomizer.getRandomDouble()-0.50) * _windowSize;
        while(newRate > 1 || newRate < 0) {
            if(newRate > 1.0) newRate = 2.0 - newRate;
            if(newRate < 0.0) newRate = -newRate;
        }
        setRate(newRate);

        _violate = false;
        return 0;
    }

    @Override
    public void undo() {
        setRate(_oldRate.doubleValue());
    }

    @Override
    public String getName() {
        return "Change-Inheritance";
    }

    @Override
    public void optimize(final double logAlpha) {
        _windowSize *= Math.exp(Utils.calcDelta(this, logAlpha));
    }

    private void setRate(double rate) {
        _target.Item1.setParentProbability(_target.Item2, rate);
        _target.Item1.setParentProbability(_target.Item3, 1.0 - rate);
    }

    // test
    public static void main(String[] args) {
        {
            UltrametricNetwork net = new UltrametricNetwork("(O:0.5,((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0:3.0);");
            SpeciesNetPriorDistribution prior = new SpeciesNetPriorDistribution(1);
            System.out.println(prior.logPrior(net.getNetwork()));
            System.out.println(net.getNetwork().toString());
            int runs = 10000;
            int test = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new ChangeInheritance(net);
                double logHR = op.propose();
                if(prior.isValid(net.getNetwork())) test++;
            }
            System.out.println(net.getNetwork().toString());
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
        }
        {
            UltrametricNetwork net = new UltrametricNetwork("(O:0.5,((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0:3.0);");
            int runs = 10000;
            int test = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new ChangeInheritance(net);
                double logHR = op.propose();
                op.undo();
                for(NetNode<NetNodeInfo> node : net.getNetwork().getNetworkNodes()) {
                    if(node.getParentProbability(node.getParents().iterator().next()) != 0.50) {
                        test--;
                        break;
                    }
                }
                test++;
            }
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
        }
    }
}
