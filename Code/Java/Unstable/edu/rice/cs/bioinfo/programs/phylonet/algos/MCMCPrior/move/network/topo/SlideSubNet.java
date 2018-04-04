package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.topo;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.distribution.SpeciesNetPriorDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by wendingqiao on 3/19/16.
 * Slide a random chosen subnet of the network
 */
public class SlideSubNet extends NetworkOperator {

    private Double _logHR;
    private NetNode<NetNodeInfo> _v1, _v2, _v3, _v4, _v5 = null, _v6 = null;
    private double _oldHeight, newHeight;
    private double _windowSize = Utils._TIME_WINDOW_SIZE;

    public SlideSubNet(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        // reset
        _logHR = null;
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null;

        List<NetNode<NetNodeInfo>> internalTreeNodes = Networks.getInternalTreeNodes(_network.getNetwork());
        _v1 = internalTreeNodes.get(Randomizer.getRandomInt(internalTreeNodes.size()));

        _oldHeight = _v1.getData().getHeight();
        newHeight = _oldHeight + getDelta();

        List<NetNode<NetNodeInfo>> children = IterableHelp.toList(_v1.getChildren());
        _v2 = children.get(Randomizer.getRandomInt(children.size()));

        _v3 = _v1.isRoot() ? null : _v1.getParents().iterator().next();
        _v4 = Networks.getOtherChild(_v1, _v2);
        double tUpperBound = _v1.isRoot() ? Double.MAX_VALUE : _v3.getData().getHeight();
        double tLowerBound = Math.max(_v2.getData().getHeight(), _v4.getData().getHeight());

        if(tLowerBound <= newHeight && newHeight <= tUpperBound) {
            // topology doesn't change, change height
            setNodeHeight(_v1, newHeight);
            _violate = newHeight > _oldHeight;
            _logHR = 0.0;
        } else if(_v4.hasParent(_v3) || newHeight <= _v2.getData().getHeight()) {
            _violate = false;
            _logHR = Utils.INVALID_MOVE;
        } else {
            Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge;
            int numDest;
            boolean downwards = newHeight < _oldHeight;

            NetNode<NetNodeInfo> root = _network.getNetwork().getRoot();
            if(newHeight >= root.getData().getHeight()) {
                numDest = 1;
                edge = new Tuple<>(root, null);
            } else {
                // downwards => search children, else search parents
                List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = downwards ?
                        getEdgesFromNodeGivenHeight(_v4, newHeight, true) :
                        getEdgesFromNodeGivenHeight(_v3, newHeight, false);
                numDest = edges.size();
                if(numDest == 0) {
                    _violate = false;
                    _logHR = Utils.INVALID_MOVE;
                    return _logHR;
                }
                edge = edges.get(Randomizer.getRandomInt(numDest));
            }
            _v5 = edge.Item2;
            _v6 = edge.Item1;

            if(_v3 != null && _v5 != null) {
                moveTail(newHeight, _v3, _v4, _v5, _v6);
            } else if (_v3 == null) {
                moveRoot(newHeight, _v4, _v5, _v6);
            } else {
                setNewRoot(newHeight, _v3, _v4, _v6);
            }
            // downwards => search parents, else search children
            int numSrc = 1;
            if (_oldHeight < _network.getNetwork().getRoot().getData().getHeight()) {
                numSrc = downwards ? getEdgesFromNodeGivenHeight(_v5, _oldHeight, false).size() :
                        getEdgesFromNodeGivenHeight(_v6, _oldHeight, true).size();
            }
            _violate = true;
            _logHR = Math.log((double) numDest / (double) numSrc);
        }

        return _logHR;
    }

    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        if(_v5 == null && _v6 == null) {
            setNodeHeight(_v1, _oldHeight);
            return;
        }
        if(_v3 != null && _v5 != null) {
            moveTail(_oldHeight, _v5, _v6, _v3, _v4);
        } else if (_v3 == null) {
            setNewRoot(_oldHeight, _v5, _v6, _v4);
        } else {
            moveRoot(_oldHeight, _v6, _v3, _v4);
        }
    }

    @Override
    public String getName() {
        return "Slide-SubNet";
    }

    private double getDelta() {
        return (Randomizer.getRandomDouble() - 0.5) * _windowSize;
    }

    @Override
    public void optimize(final double logAlpha) {
        _windowSize *= Math.exp(Utils.calcDelta(this, logAlpha));
    }

    private void setNewRoot(double height,
                            NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                            NetNode<NetNodeInfo> v6) {

        _v1.getData().setHeight(height);
        _v2.setParentDistance(_v1, height - _v2.getData().getHeight());

        double rootPopSize = v6.getRootPopSize();
        double[] paramV3V1 = getParameters(v3, _v1);
        double[] paramV1V4 = getParameters(_v1, v4);

        v3.removeChild(_v1);
        _v1.removeChild(v4);

        adopt(v3, v4, paramV1V4);
        adopt(_v1, v6, paramV3V1);

        _network.getNetwork().resetRoot(_v1);
        _v1.setRootPopSize(rootPopSize);
    }



    private void moveRoot(double height,
                          NetNode<NetNodeInfo> v4,
                          NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6) {

        _v1.getData().setHeight(height);
        _v2.setParentDistance(_v1, height - _v2.getData().getHeight());

        double rootPopSize = _v1.getRootPopSize();
        double[] paramV1V4 = getParameters(_v1, v4);
        double[] paramV5V6 = getParameters(v5, v6);

        _v1.removeChild(v4);
        v5.removeChild(v6);

        adopt(v5, _v1, paramV1V4);
        adopt(_v1, v6, paramV5V6);

        _network.getNetwork().resetRoot(v4);
        v4.setRootPopSize(rootPopSize);
    }

    private void moveTail(double height,
                          NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                          NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6) {

        _v1.getData().setHeight(height);
        _v2.setParentDistance(_v1, height - _v2.getData().getHeight());

        double[] paramV3V1 = getParameters(v3, _v1);
        double[] paramV1V4 = getParameters(_v1, v4);
        double[] paramV5V6 = getParameters(v5, v6);

        v3.removeChild(_v1);
        _v1.removeChild(v4);
        v5.removeChild(v6);

        adopt(v3, v4, paramV1V4);
        adopt(v5, _v1, paramV3V1);
        adopt(_v1, v6, paramV5V6);
    }

    private List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>>
    getEdgesFromNodeGivenHeight(NetNode<NetNodeInfo> node, double height, boolean downwards) {

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> list = new ArrayList<>();

        Stack<NetNode<NetNodeInfo>> stack = new Stack<>();
        stack.add(node);
        Set<NetNode<NetNodeInfo>> visited = new HashSet<>();

        while(!stack.isEmpty()) {
            NetNode<NetNodeInfo> n = stack.pop();
            if(visited.contains(n)) continue;
            visited.add(n);
            if(downwards) { // search children
                for(NetNode<NetNodeInfo> child : n.getChildren()) {
                    if(child.getData().getHeight() < height) {
                        if(child != _v2) {
                            list.add(new Tuple<>(child, n));
                        }
                    } else {
                        stack.add(child);
                    }
                }
            } else { // search parent
                for(NetNode<NetNodeInfo> par : n.getParents()) {
                    if(par.getData().getHeight() > height) {
                        list.add(new Tuple<>(n, par));
                    } else {
                        stack.add(par);
                    }
                }
            }
        }
        return list;
    }

    public String printMove() {
        return _v1.getName() + _v2.getName()
                + (_v3 == null ? "null" : _v3.getName()) + _v4.getName()
                + (_v5 == null ? "null" : _v5.getName()) + _v6.getName()
                + "  " + _oldHeight + "  " + newHeight;
    }

}
