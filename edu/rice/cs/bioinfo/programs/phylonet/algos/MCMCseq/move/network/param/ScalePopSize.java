package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.param;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by wendingqiao on 3/3/16.
 * Multiply a random population size by a random factor.
 */
public class ScalePopSize extends NetworkOperator {

    private double _scaleFactor = 0.92;
    private double _upperLimit = 1.0 - 1e-6;
    private double _lowerLimit = 0.8;

    private Double scale;
    private List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> _edges;
    private List<Double> _prevPopSizes;
    private double _prevRootPopSize = -1;

    public ScalePopSize(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _edges = null;
        scale = null;

        scale = getScaler();

        int dimension = 1;
        if(Utils._CONST_POP_SIZE) {
            _edges = null;
            _prevPopSizes = null;
        } else {
            _edges = Networks.getAllEdges(_network.getNetwork());
            dimension += _edges.size();
        }
        scalePopSize(scale, false);

        _violate = false;
        return Math.log(scale) * (dimension - 2);
    }

    @Override
    public void undo() {
        scalePopSize(1.0 / scale, true);
    }

    @Override
    public String getName() {
        return "Scale-PopSize";
    }

    @Override
    public void optimize(double logAlpha) {
        double delta = Utils.calcDelta(this, logAlpha);
        delta += Math.log(1.0 / _scaleFactor - 1.0);
        setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
    }

    private void setCoercableParameterValue(final double value) {
        _scaleFactor = Math.max(Math.min(value, _upperLimit), _lowerLimit);
    }

    private double getScaler() {
        return _scaleFactor + Randomizer.getRandomDouble() * (1.0 / _scaleFactor - _scaleFactor);
    }

    private void scalePopSize(double scale, boolean undo) {
        NetNode root = _network.getNetwork().getRoot();
        if(undo) {
            root.setRootPopSize(_prevRootPopSize);
            _prevRootPopSize = -1;
        } else {
            _prevRootPopSize = root.getRootPopSize();
            root.setRootPopSize(_prevRootPopSize * scale);
        }
        if(Utils._CONST_POP_SIZE) {
            return;
        }
        if(_edges == null || (undo && _prevPopSizes == null)) {
            throw new RuntimeException("not CONST population size!!!!");
        }
        if(undo) {
            for(int i = 0; i < _edges.size(); i++) {
                _edges.get(i).Item1.setParentSupport(_edges.get(i).Item2, _prevPopSizes.get(i));
            }
            _prevPopSizes = null;
        } else {
            _prevPopSizes = new ArrayList<>();
            for(Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge : _edges) {
                double ps = edge.Item1.getParentSupport(edge.Item2);
                _prevPopSizes.add(ps);
                edge.Item1.setParentSupport(edge.Item2, ps * scale);
            }
        }
    }

    public static void main(String[] args) {
        Utils._CONST_POP_SIZE = false;
        {
            UltrametricNetwork net = new UltrametricNetwork("(((AF:4.7560289410455474E-7:6.216710846747218E-4)I5#H1:6.673080055619683E-5:6.04717164578389E-4:0.19169217458789256,EU:6.720640345030139E-5:4.54804059863962E-4)I4:3.2385691593548525E-5:1.441308148785055E-4,(I5#H1:7.170661624805849E-5:8.149685305377055E-4:0.8083078254121074,SA:7.218221914216304E-5:0.0011303870974394095)I1:2.740987590168687E-5:3.0250110970208124E-4)I0;");
            System.out.println(net.getNetwork().toString());
            int runs = 10000;
            int counter = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new ScalePopSize(net);
                double logHR = op.propose();
                if(logHR != Utils.INVALID_MOVE) counter++;
            }
            System.out.println(counter == runs);
            System.out.printf("%d out of %d\n", counter, runs);
            System.out.println(net.getNetwork().toString());
        }
        {
            UltrametricNetwork net = new UltrametricNetwork("(((AF:4.7560289410455474E-7:6.216710846747218E-4)I5#H1:6.673080055619683E-5:6.04717164578389E-4:0.19169217458789256,EU:6.720640345030139E-5:4.54804059863962E-4)I4:3.2385691593548525E-5:1.441308148785055E-4,(I5#H1:7.170661624805849E-5:8.149685305377055E-4:0.8083078254121074,SA:7.218221914216304E-5:0.0011303870974394095)I1:2.740987590168687E-5:3.0250110970208124E-4)I0;");
            System.out.println(net.getNetwork().toString());
            int runs = 10000;
            int test = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new ScalePopSize(net);
                double logHR = op.propose();
                op.undo();
                if(net.getNetwork().getRoot().getRootPopSize() != Utils._POP_SIZE_MEAN) {
                    test--;
                }
                test++;
            }
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
            System.out.println(net.getNetwork().toString());
        }
    }
}
