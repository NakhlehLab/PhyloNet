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
 * Created by wendingqiao on 3/7/16.
 * select a random internal node, select a random child of node
 * select a random height between [child, node]
 * select a random edge cross with the height where edge.child != child and edge.parent != node
 * swap parents
 * the same as SMARTIE
 * http://sysbio.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=20525618
 */
public class SwapNodes extends NetworkOperator {

    private double _logHR;
    private Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> _edge1, _edge2;

    public SwapNodes(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _edge1 = _edge2 = null; //reset

        List<NetNode<NetNodeInfo>> internalNodes = Networks.getInternalNodes(_network.getNetwork());
        NetNode<NetNodeInfo> parent, child;
        do {
            parent = internalNodes.get(Randomizer.getRandomInt(internalNodes.size()));
        } while (parent.isRoot());

        List<NetNode<NetNodeInfo>> children = IterableHelp.toList(parent.getChildren());
        child = children.get(Randomizer.getRandomInt(children.size()));
        _edge1 = new Tuple<>(child, parent);

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = getEdgesGivenEdgeHeight(_edge1);
        if(edges.size() == 0) {
            _logHR = Utils.INVALID_MOVE;
            _violate = false;
        } else {
            _edge2 = edges.get(Randomizer.getRandomInt(edges.size()));

            if(_edge1.Item1.hasParent(_edge2.Item2) || _edge2.Item1.hasParent(_edge1.Item2)) {
                _logHR = Utils.INVALID_MOVE;
                _violate = false;
            } else {
                swap(_edge1.Item1, _edge1.Item2, _edge2.Item1, _edge2.Item2);

                _logHR = 0;
                _violate = true;
            }
        }
        return _logHR;
    }

    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        swap(_edge1.Item1, _edge2.Item2, _edge2.Item1, _edge1.Item2);
    }

    @Override
    public String getName() {
        return "Swap-Nodes";
    }

    private void swap(NetNode<NetNodeInfo> c1, NetNode<NetNodeInfo> p1, NetNode<NetNodeInfo> c2, NetNode<NetNodeInfo> p2) {

        double rootPopSize = _network.getNetwork().getRoot().getRootPopSize();
        if(rootPopSize == Double.NaN) throw new RuntimeException("Invalid root pop size in SwapKernel!!!");

        double[] paramP1C1 = getParameters(p1, c1);
        double[] paramP2C2 = getParameters(p2, c2);

        p1.removeChild(c1);
        p2.removeChild(c2);

        adopt(p1, c2, paramP2C2);
        adopt(p2, c1, paramP1C1);
    }

    private List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>>
    getEdgesGivenEdgeHeight(Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge) {

        double height = edge.Item2.getData().getHeight();
        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = new ArrayList<>();
        for(NetNode<NetNodeInfo> child : _network.getNetwork().dfs()) {
            if(child.getData().getHeight() >= height || child.hasParent(edge.Item2)) continue;
            for(NetNode<NetNodeInfo> par : child.getParents()) {
                if(par.getData().getHeight() < height) continue;
                edges.add(new Tuple<>(child, par));
            }
        }
        return edges;
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
            int test = 0;
            int test1 = 0;
            int test2 = 0;
            for(int i = 0; i < runs; i++) {
                SwapNodes op = new SwapNodes(net);
                double logHR = op.propose();
                op.undo();
                test++;
                if(Networks.hasTheSameTopology(topo, net.getNetwork())) {
                    test1++;
                } else {
                    System.out.println(net.getNetwork().toString());
                    System.out.println(op._edge1.Item1.getName() + " " + op._edge1.Item2.getName());
                    System.out.println(op._edge2.Item1.getName() + " " + op._edge2.Item2.getName());
                    System.exit(1);
                }
                if(prior.isValid(net.getNetwork()) && net.isUltrametric()) test2++;
            }
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
            System.out.println(test1 == runs);
            System.out.printf("%d out of %d\n", test1, runs);
            System.out.println(test2 == runs);
            System.out.printf("%d out of %d\n", test2, runs);
            System.out.println(net.getNetwork());
        }
        {
            Map<String, Double> heights = new HashMap<>();
            for(NetNode<NetNodeInfo> node : Networks.postTraversal(net.getNetwork())) {
                heights.put(node.getName(), node.getData().getHeight());
            }
            int runs = 10000;
            int counter0 = 0, counter1 = 0, counter2 = 0, counter3 = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new SwapNodes(net);
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
            for(NetNode<NetNodeInfo> node : Networks.postTraversal(net.getNetwork())) {
                if(Math.abs(heights.get(node.getName()) - node.getData().getHeight()) > 0.000001) {
                    System.out.println(node.getName() + " before " + heights.get(node.getName()) + " after " + node.getData().getHeight());
                    break;
                }
            }
            System.out.printf("%d out of %d valid move!!!\n", counter1, runs);
            System.out.println(counter2 == runs);
            System.out.printf("%d out of %d valid network!!!\n", counter2, runs);
        }
    }

}
