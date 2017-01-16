package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.topo;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SpeciesNetPriorDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.start.UPGMATree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.start.distance.JCDistance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
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

    public static void main(String[] args) {
        Map<String, String> locus = new HashMap<>();
        {
            locus.put("A", "CTTCGTGACGGGCTCGGCTCGTACGGTCAAGGGCACCTGAGCTAGGCAACTCAAGACGGGCGAGAGTCCCGCTACTACGGAAAAGGGTTTGCAGCTTCGAACATAGAACCACGGGACCTCCGGGACTCGCCCAATTGCGGTCGATGGCACTACGATTGGACAGGCGCTTACGCCAATTTACAGGTAGCGAGAGGTCTGCGCAGCAAAAGACTCCGCTTCACCCGGTGGTACCTGACTCGCGGCGCTGACCGTTCGCCGTAAGCGTTGACAATTTTCGGAAGCATCGCCAAAGTGCTGGAGTACGGGCTTCACACCCGCAACCACCGACAGCCCAACGAATCACTTCGCCGGTTGTCCTTCACCCCTGTTGGCAGAGGATGCTTGCGTGACTTTCATCCTTCTGTCCTTTCGTGCGGACTCGCACAGATCCTCCAAACAAGCGAGATCCGACCGATACTCTGCCCCCAGCAAGGCGGGGTTTCAGCGTCCCGTAGCTAGTAGGTATCCAGTCCAGAGGCACGTAGAGCACTCACCGTCCCCAGCCCCCCTCCACTCATCGCTGTGTTGAGGTCAAGTTCCCCGAGATCGTACTAACCGGTATGCAACCCACATCTGCCCATAGTCTACCCATCATACTACAGTCAGTAAACCCGGAGTGTATGGCTCGACTAATGTCCGTACAGCCAGCGCTCGATAAGTGCGTGCTCGAACGCTTACACCCCTGCTGGCACCAGGTAGCGCATGATCCACCAATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCGAAAGTATCGAGTTAAGCCTCAGCAATATGGGCCACCTTGACTGAGCTACACTCCCTCCCGTACGGGTGAAGTCCCCGCCGGCACAGGGGGGACACTTCTATTGAACATTTCGCCTCACCGAAGTCGCCCCCGAGTACGCCAGTGACGTACAGTCGCACGCGGGTTTGGTAAGGAGGAGCCAGATCGCGTCTACATGTTGAGTAGGTCCCGGCGC");
            locus.put("C", "CTACGTGACGGGCCCGGCTCGCACGGTCAATGGCACCTGCGCTAGGCAACTCGAGACGGGCGAGGGTCCCGCTACTACGGACAAGGGTGTGCAGCGTCGAACATAGAACCACGGGACCTCCGGTACTCGCTCAACCGCGGTCGGTGGTACTACGATTGGACGGGCACTTTCGGCAATTTACAGTTAGCGAGAGGTCTGCGCGGCAAAAGACTCCAATTCACCCGGCGTTACCTGACTCGCGGCACTGACCGTTCGCCGTAAGCGTTGACAATTTTCGGTAGCATCGCCAAAGTGCTGGAGTACGGGTTTCACACTTGTAGCCACCGACAGCCCAACGGATCAATTCGCCAGTTGCCCTTCACTCCTGTTAGCAGAGGATGCTAGCGTGAGTTTCATCCATCTGTCCTGTCGAGCGGACTCGCACAGATACTCCAAACAAGTGAGATCCGTCCGGTACTCTACCCCTCACAAAGCGGGATTTCAGCGTCTCGTAGTGAGTAGGCATCCGGTCCAGGGGCACGTAGAGCACTCACCGTCCTCAGCCCCCTTCTATTCATTGCTGTGTTGAAGCGAAGTTCCCCGAAATCGTACTAACCGGTATGCAACCCAGATCTGCCCATAGTTTCCACATCATACTACAGCCAGAAAACTCGGAGTGTATGGCTTGACTAATGTCCGTACTGCCAGCGCTAGGTAAGTACGTGCTCGAACGCTTACACCCCTGCTGGCACCTGGTTGCGCATGATCCACCAATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCAAAAGTATTGAATTAAGCCTCAGCAATATGGGCCACCCAGACTGAGCCACGCTCCCCTCCGTACGGATGAAGTCCCCACCGGGGCAGGGGGGACGCTTCTATTGAACACTTCGCCTCACCGAAGTCGCCCCCGAGTACGCCAGTGACGTACGGTCGCACGCGGGTTCGGTAAGGAGGGGTCAGATCGCGTCTGCATGTTGAGTAGGTCCCTGCGC");
            locus.put("G", "CTACGTGACGGGCCTGGCTCGCACGGTCAATGGCACTTGAGCTAGGCAACCCAAGACGGGCGAGGGTCCCGCTACTACGGAAAAGGGTGTGCAGCTTCGAACATAGAACCACGGGACCTCCGGGACTCGCTCAATTGCGGTCGGTGGTACTACGATTGGACGGGCACTTTCGGCAATTTACAATTAGCGAGAGGTCTGCGCAGCAAAAGACTCCAATTCACCCGGCGCTACCTGATTCGCGGCACTGACCGTTCGCCGTAAGCGTTGACAAATTTCGGTAGCATCGCCAAAGTGCTGGAGTACGGGCTTCACACTTGTAGCCACCGACAGCCCAACGAATCAATTCGCCAGTTGCCCTTCACTCCTGTTAGCAGAGGATGCCAGCATGAGATTCATCCATCTGTCCTTTCGAGCGGACTCGCACAGATACTCCAAACAAGCGAGATCCGACCGGTACTCTACCCCTCACAAAGCGGGATTTCAGCGTCTCGTAGTGAGTAGGTATCCGGTCCAGAGGCACGTAGAGCCCTCACTGTCCTCAGCCCCCCTCTACTCATTGCTGTGTTGAGGTGAAGTTCCCCGAAATCGTACTAACCGGTATGCAACCCAGATCTGCCCATAGGCTACACATCATACTACAGCCAGTAAACCCGGAGTGTATGGCTTGACTAAGGTCCGTACAGCCAGCGCTCGATAAGTACGTGCTCGAACGCTTACACCCCTGCTGGCACCAGGTAGCGTATGATCCACCGATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCGAAAGTATTGAGTTAAGCCTCAGCAATATGGGCCACCCTGACTGAGCTACGCTCCCTTCCGTACGGATGAAGTCCCCACCGGGACAGGGGGGACGCTTCTATTGAACACTTCACCTCACCGAAGTCGCCCCCGAGTACACCAGTGACGTACAGTCGCACGCGGGTTCGATAAGGAGGGGCCAGATCGCGTCTACATGTTGAGTAGGTCCCTGCGC");
            locus.put("R", "TTGTGCGACGGGCGCGGCCCGGGCGATCAAGGGTACCAGAGTCAGGCAACTTAAGATGGGCGAGAGCCCCGCTAATACGGAAAGGAGTATGCAGCTCCGGACATGGGATCACTGGATCTCCAGGACTCGGTCGACTGCGGTCGGCGTTACTACAATCGGAGAGGCACATTCGCTAACTCATAGTTTGCGAGAGATCTGCGGAGCAAAGGATTCCCCTCTACTCGGCGCTACCTGACTCGCTACGCAGACCGTTCGCCGTAAGTGTTGCCAATCCCCGGAGGCATCGCCAAAGTACTGGAGTACGGGCTTAACACCTACAGCCACCGACAGCGCAACGAATCATTTCACCAGTTGCCTTTTACCCTTGTTAGCTGAGGATGCTAGCTTAACTTTCATCCTATGTTCCTCTCGGGCGGACTCTAACTGATCCCCCGAATGAGCGAGGTCTGACCGATACTCTGCCCCAGGCAAGGGGGGAGCCCATCCCCTTATAGTAAGCAAATCCCCAGTTCAGAGGCACATAGAGCACCCACCGTCCACAGCCCCTTTCCACTCACTGGTGCGCTGAGGTGAAAGTGCCCGAAATCCTACGAATTGGTATGCAACCCAGATCTGTAGGCAGGCTACATGTCATACTACAGCTCGTAAACTTGGAGTGTATGGCTAGACTGATATCCGAACAAACAACGCTCGACAAGCGCGACCTCGACCGCTCACACCCTTGCTGACACCAAGCAACACATGATCCATCAGTGCAGCCCCAACGTTTTTTGTGACCTCCGTCCGAAAGTATGGATTTGAGCCTCAGCAATGTGGCCCACCATGGCCGAGCTACGCTCCCCTACGTACGGATGATTTCCCCGCCGGGACAGGCGGGACGGTTCTATTAAACATCTCACCTTACTGATGTCGCCCCCGGGTACGGCAGCGACGTACAGCCGCACGCGAGCTTGGTAAGGAGGAGCCAGATCGTGCCTACATGTTGAGTAGGTCCCTACTC");
            locus.put("Q", "CCATACGATGGGCTCGGCTCGTATGATTGAGGGCACCGAAGCTAGGCGACTCAAGATGGGCGAGGGCTCCGCGAATACGGAAAAGGGTATGCAGCTTCGGCCATAGGACCACGTGATCTCCGGGACTCGCCCAATTGAGATCGGCGTTACTACAATTAGACGAGCACATTCGCCAATTTATAGTTAACAAGAGATCTGCGTAGCAAAAGATTCAACGTTACCCGGCGCTACCAGACCCGCGGCACAGGCCGTTTGCCCTAAGCGTTGACAATCTTCGGATTCATCGCCAAAGTGCTGGAGTACAGGCTTCACACCTACAGCCACTGACAGCCTAACGAATCACTTCACCAATTGCCTTTCACCCCTGTTAGCGGAGGGCGCTAGCATAACTTTCGTCCTACGTTCCTCTCGTACGGATTCGGACAGATCCTCCGAGCAAGCGAAGTCCGGCCGATACTCTGCCCCTAGCAAGGCGGGATCCCGTCGCCTTGTAGTGAGCAGATATCCAGTTTGGGAGCACATAGAGCACCGACCGTCCACAGTCCCCTTCTCTTCATTGGTGCGTTGAGGTGAAAGTTCCCGAAATCCTACAAACTGGTATGCAAACCGGATCTGCAGGTTGGCTACGTATCATACTACAGCCCGTAAACTCGGAGTGCATGGCTTGACTAACATCCGTACAAACAGCGCTCGATGAGTGCGACCTCGGCAGCTTACACCCTTGCTGACACCAAATAGCGCATGATCCACCAGTACAGCCCAAGCGCCTTTCGCGTCCTCCGCCCGAATGTATGGATGTAAGCCTGAGTCACGTGGACCACCGTGCCCGAGCTACGCTCCCTTACGTGCGGATGATGTCCCCGCCGGGACAGGCGGAACGCTTCTATTGAACATTTCACCTCGCTGAAGTCGCCCCTGAGTACGCCAGCAACGTACAGTCGCACGCGAGCTCAGTAAGGAGGAGCCATATCGCGTCTACATGATGCGTAGGTCCCTGCGC");
            locus.put("L", "CTACGCGACGGGCTCGGCTCGTACAATAAAGGGCACCAGAGCTAGGTAACCCAAGATAGGCGAGGGCTTCGCTAATACGGAGAAGGGTATTCAGCTTCGGTCACAGGGCCACGTGATCTCCGGGACTCGCCCAGTTGCGATCGGCGTTACTACAATTGGAGGGGCACATTCGCCAACTTACAGTTAGTGAGGGATCTGCGTATCAAAAAACTCCACACTACTCGGCGCGACCTGGCTCGTGGCGCAGGCCGTTTGCCGTAAGTGTTGAGAATTTTCGGAAGAATGGCCGGAGTGTTAGAGTCCAGGCTTCACACCTACAGCCATTGACAGCTCAACGAGTCATTTCACCAGTTGCCTTTCACCCCTGTTAGCAGAGGGAGCTAGCGCAGCTTTCATCCTATGTCCATCTCGGGCGGACTCGGACAGATCCTCCGAACAAGCGAGGTCCGGTCGATACTCTGCCCCTAACAAGGCAGGATCCCATCGCCTTGTAGTGAGCAGATATCCAGTTTAGGGGCGCATAGAGCACCCACCGTCCACAGCCCCCTTGTACTTATTGGTGTGTTGAGGTGAAAGTCCCCGAAATCCCAGGAACTTGTATGCAACCCAGATCTGCAGGTAGGCTACGTATCATACTGTACCCCGCAAACTCGGAGTGTGTGACTGGACTAACGTCCGTACAAACAGCGCTCGATAGGTGCGACCTCGACAGCTTACACCCTTGCTGGCACCAAATAGCGCGTGATCCACCAGTGCAGACCAAACGTTTTCTGCGCCCTCCGTCCGAAAGTACGGGTTTAAGTCTCAGTAACATGGTCCACTATGACCGAACTATACCCCCTTACGTACGGGTGATGTCCCCGCCGGGACGGGCGGGACGCTTCTACTGAACGCTTCACTTCACTGAAGTCGCCCCTGAGTATGCCAGCAACGTACAGTCGCACGCGAGCTCGGTAAAGAGGAACCAGATCCCGTCGACATGTTGAGTTGGTCCCTATGC");
        }
        Alignment seq = new Alignment(locus);
        UPGMATree upgma = new UPGMATree(new JCDistance( seq.getAlignment() ));
        List<UltrametricTree> geneTrees = new ArrayList<>();
        geneTrees.add(upgma.getUltrametricTree());

        SpeciesNetPriorDistribution prior = new SpeciesNetPriorDistribution(1);
        UltrametricNetwork net = new UltrametricNetwork("((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0:3.0;"
                , geneTrees);
        System.out.println(net.getNetwork());
        System.out.println(prior.logPrior(net.getNetwork()));
        Network topo = Networks.readNetwork(net.getNetwork().toString());
        {
            int runs = 10000;
            int counter = 0;
            int test1 = 0;
            int test2 = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new SlideSubNet(net);
                if(op.propose() != Utils.INVALID_MOVE) {
                    counter++;
                }
                prior.isValid(net.getNetwork());
                op.undo();
                if(Networks.hasTheSameTopology(topo, net.getNetwork())) {
                    test1++;
                } else {
                    System.out.println("after undo - " + net.getNetwork().toString());
                    System.exit(1);
                }
                if(prior.isValid(net.getNetwork()) && net.isUltrametric()) {
                    test2++;
                } else {
                    System.err.println(net.getNetwork().toString());
                    System.exit(1);
                }
            }
            System.out.printf("%d out of %d\n", counter, runs);
            System.out.printf("%d out of %d\n", test1, runs);
            System.out.println(test1 == runs);
            System.out.printf("%d out of %d\n", test2, runs);
            System.out.println(test2 == runs);
            System.out.println(net.getNetwork());
        }
        {
            int runs = 10000;
            int counter0 = 0, counter1 = 0, counter2 = 0, counter3 = 0;
            for(int i = 0; i < runs; i++) {
                SlideSubNet op = new SlideSubNet(net);
                Network beforeProposal = Networks.readNetwork(net.getNetwork().toString());
                if(!net.isValid()) {
                    System.err.println("invalid temporal constraints before proposal!!!!");
                    System.exit(1);
                }
                double logHR = op.propose();
                String afterProposal = "empty";
                if(prior.isValid(net.getNetwork())) {
                    afterProposal = net.getNetwork().toString();
                    if(net.isValid() && logHR != Utils.INVALID_MOVE) {
                        if(net.isUltrametric()) counter0++;
                        counter1++;
                    } else {
                        op.undo();
                    }
                    counter2++;
                } else {
                    op.undo();
                }
                if(!prior.isValid(net.getNetwork()) || !net.isUltrametric() || !net.isValid()) {
                    System.err.println("!!!!!!!!! wrong Invaild move after undo \n" + beforeProposal.toString() + "\n"
                            + afterProposal + "\n"
                            + Networks.hasTheSameTopology(beforeProposal, net.getNetwork()));
                    System.err.println(net.getNetwork().toString());
                    System.err.println(op.printMove());
                    System.exit(1);
                }
                if(Networks.hasTheSameTopology(beforeProposal, net.getNetwork())) {
                    counter3++;
                }
            }
            System.out.println(net.getNetwork());
            System.out.println(counter0 == counter1);
            System.out.printf("%d + %d = %d ~ %d\n", counter1, counter3, counter1 + counter3, runs);
            System.out.println(!Networks.hasTheSameTopology(net.getNetwork(), topo));

            for(NetNode<NetNodeInfo> node : net.getNetwork().bfs()) {
                for(NetNode<NetNodeInfo> par: node.getParents()) {
                    if(node.getParentSupport(par) != Utils._POP_SIZE_MEAN) {
                        System.out.println("popSize: " + node.getName() + " before " + Utils._POP_SIZE_MEAN
                                + " after " + node.getParentSupport(par) );
                    }
                }
            }
            System.out.printf("%d out of %d valid move!!!\n", counter1, runs);
            System.out.println(counter2 == runs);
            System.out.printf("%d out of %d valid network!!!\n", counter2, runs);
        }
    }

}
