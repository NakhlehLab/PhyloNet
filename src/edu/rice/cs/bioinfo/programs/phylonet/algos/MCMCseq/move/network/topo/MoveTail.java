package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.topo;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution.SpeciesNetPriorDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.UPGMATree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.distance.JCDistance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 3/16/16.
 * Move the tail of a random selected edge to a random source
 */
public class MoveTail extends NetworkOperator {

    private double _logHR;
    private NetNode<NetNodeInfo> _v1, _v2, _v3, _v4, _v5, _v6;
    private double _oldHeight;

    public MoveTail(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null; // reset

        List<NetNode<NetNodeInfo>> internalTreeNodes = Networks.getInternalTreeNodes(_network.getNetwork());
        do {
            _v1 = internalTreeNodes.get(Randomizer.getRandomInt(internalTreeNodes.size()));
        } while(_v1.isRoot());

        List<NetNode<NetNodeInfo>> children = IterableHelp.toList(_v1.getChildren());
        _v2 = children.get(Randomizer.getRandomInt(children.size()));

        _v3 = _v1.getParents().iterator().next();
        _v4 = Networks.getOtherChild(_v1, _v2);

        if(_v4.hasParent(_v3)) {
            _violate = false;
            _logHR = Utils.INVALID_MOVE;
        } else {
            List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = getEdgesAboveNode(_v2, _v1);
            if(edges.size() == 0) {
                _violate = false;
                _logHR = Utils.INVALID_MOVE;
            } else {
                _oldHeight = _v1.getData().getHeight();

                Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge = edges.get(Randomizer.getRandomInt(edges.size()));
                _v5 = edge.Item2;
                _v6 = edge.Item1;

                double t5 = _v5.getData().getHeight();
                double tl = Math.max(_v2.getData().getHeight(), _v6.getData().getHeight());
                double newHeight = tl + Randomizer.getRandomDouble() * (t5 - tl);

                moveTail(newHeight, _v3, _v4, _v5, _v6);

                double l2 = t5 - tl;
                double l1 = _v3.getData().getHeight() - Math.max(_v2.getData().getHeight(), _v4.getData().getHeight());

                _violate = true;
                _logHR = Math.log(l2 / l1);
            }
        }
        return _logHR;
    }

    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        moveTail(_oldHeight, _v5, _v6, _v3, _v4);
    }

    @Override
    public String getName() {
        return "Move-Tail";
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

    private List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> getEdgesAboveNode( NetNode<NetNodeInfo> v2
            , NetNode<NetNodeInfo> v1 ) {

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> list = new ArrayList<>();
        for(NetNode<NetNodeInfo> node : _network.getNetwork().bfs()) {
            for(NetNode<NetNodeInfo> par: node.getParents()) {
                if(par.getData().getHeight() > v2.getData().getHeight()
                        && par != v1 && node != v1 && node != v2) {
                    list.add(new Tuple<>(node, par));
                }
            }
        }
        return list;
    }

    // test
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
                NetworkOperator op = new MoveTail(net);
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
                NetworkOperator op = new MoveTail(net);
                Network beforePropose = Networks.readNetwork(net.getNetwork().toString());
                double logHR = op.propose();
                if(prior.isValid(net.getNetwork())) {
                    if(net.isValid() && logHR != Utils.INVALID_MOVE) {
                        if(net.isUltrametric()) counter0++;
                        counter1++;
                    } else {
                        op.undo();
                    }
                    counter2++;
                } else {
                    System.out.println("Invaild network move " + beforePropose);
                    System.out.println(net.getNetwork().toString());
                    op.undo();
                }
                if(Networks.hasTheSameTopology(beforePropose, net.getNetwork())) {
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
