package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.dimension;

import edu.rice.cs.bioinfo.library.programming.Tuple;
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 3/22/16.
 * Abstract operator that changes dimension of network
 */
public abstract class DimensionChange extends NetworkOperator {

    public DimensionChange(UltrametricNetwork net) {
        super(net);
    }

    protected double removeReticulation(NetNode<NetNodeInfo> v1, NetNode<NetNodeInfo> v2,
                                        NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                                        NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6) {

        double[] paramV1V4 = getParameters(v1, v4);
        double[] paramV2V6 = getParameters(v2, v6);
        double a = Math.min(paramV1V4[1], paramV2V6[1]);
        double b = Math.max(paramV1V4[1], paramV2V6[1]);
        double logProb = Utils.varyPopSizeAcrossBranches() ?
                getLogPopSizeProb(v1.getParentSupport(v3), a, b) +
                getLogPopSizeProb(v2.getParentSupport(v1), a, b) +
                getLogPopSizeProb(v2.getParentSupport(v5), a, b) : 0;

        v3.removeChild(v1);
        v1.removeChild(v4);
        v5.removeChild(v2);
        v2.removeChild(v6);
        v1.removeChild(v2);

        adopt(v3, v4, paramV1V4);
        adopt(v5, v6, paramV2V6);
        return logProb;
    }

    protected double addReticulation(NetNode<NetNodeInfo> v1, NetNode<NetNodeInfo> v2,
                                     NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                                     NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6,
                                     double gamma) {

        double[] paramV3V4 = getParameters(v3, v4);
        double[] paramV5V6 = getParameters(v5, v6);

        double a = Math.min(paramV3V4[1], paramV5V6[1]);
        double b = Math.max(paramV3V4[1], paramV5V6[1]);

        Tuple<Double, Double> popSizeV3V1 = sample(a, b);
        Tuple<Double, Double> popSizeV5V2 = sample(a, b);
        Tuple<Double, Double> popSizeV1V2 = sample(a, b);

        v3.removeChild(v4);
        v5.removeChild(v6);

        adopt(v3, v1, new double[] {NetNode.NO_PROBABILITY, popSizeV3V1.Item1} );
        adopt(v1, v4, paramV3V4);
        adopt(v5, v2, new double[] {1.0-gamma, popSizeV5V2.Item1} );
        adopt(v2, v6, paramV5V6);

        adopt(v1, v2, new double[] {gamma, popSizeV1V2.Item1} );
        return Math.log(popSizeV3V1.Item2) + Math.log(popSizeV5V2.Item2) + Math.log(popSizeV1V2.Item2);
    }

    protected void undoDeleteReticulation(NetNode<NetNodeInfo> v1, NetNode<NetNodeInfo> v2,
                                     NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                                     NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6,
                                     double gamma, double popSizeV3V1, double popSizeV5V2, double popSizeV1V2) {

        double[] paramV3V4 = getParameters(v3, v4);
        double[] paramV5V6 = getParameters(v5, v6);

        double a = Math.min(paramV3V4[1], paramV5V6[1]);
        double b = Math.max(paramV3V4[1], paramV5V6[1]);

        v3.removeChild(v4);
        v5.removeChild(v6);

        adopt(v3, v1, new double[] {NetNode.NO_PROBABILITY, popSizeV3V1} );
        adopt(v1, v4, paramV3V4);
        adopt(v5, v2, new double[] {1.0-gamma, popSizeV5V2} );
        adopt(v2, v6, paramV5V6);

        adopt(v1, v2, new double[] {gamma, popSizeV1V2} );
    }

    // generate sample & partial derivative from sampling distribution
    private Tuple<Double, Double> sample(double a, double b) {
        boolean setPopSize = Utils.varyPopSizeAcrossBranches();
        if(!setPopSize) {
            return new Tuple<>(Double.NEGATIVE_INFINITY, 1.0);
        }
        double c = 1.5 / (2 * b - a);
        double u = Randomizer.getRandomDouble();
        if(u < c * a / 2.0) {
            return new Tuple<>(Math.sqrt(2.0 * a * u / c), Math.sqrt(a / 2.0 / c / u));
        } else if (u <= 0.75) {
            return new Tuple<>(u / c + 0.5 * a, 1.0 / c);
        } else {
            return new Tuple<>(-0.25 / c * Math.log(4.0 * (1.0 - u)) + b, 0.25 / c / (1.0 - u));
        }
    }

    private double getLogPopSizeProb(double x, double a, double b) {
        double c = 1.5 / (2 * b - a);
        if(x < a) {
            return Math.log(c / a * x);
        } else if(x <= b) {
            return Math.log(c);
        } else {
            return Math.log(c) - 4.0 * c * (x - b);
        }
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
            int counter = 0;
            int test1 = 0;
            int test2 = 0;
            for(int i = 0; i < runs; i++) {
                AddReticulation op = new AddReticulation(net);
                String networkBefore = net.getNetwork().toString();
                double logHR = op.propose();
                if(prior.isValid(net.getNetwork()) && net.isUltrametric()) {
                    counter++;
                } else {
                    System.err.println("propose failure");
                }
                op.undo();
                test++;
                if(Networks.hasTheSameTopology(topo, net.getNetwork())) {
                    test1++;
                } else {
                    System.out.println(net.getNetwork().toString());
                    System.exit(1);
                }
                if(prior.isValid(net.getNetwork()) && net.isUltrametric()) {
                    test2++;
                } else {
                    System.err.println("undo failure");
                    System.err.println(networkBefore);
                    System.err.println(net.getNetwork().toString());
                    System.err.println(op._v1.getName() + " : " + op._v2.getName() + " : " + op._v3.getName() + " : "
                            + op._v4.getName() + " : " + op._v5.getName() + " : " + op._v6.getName());
                    System.exit(1);
                }
            }
            System.out.println(counter == runs);
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
            System.out.println(test1 == runs);
            System.out.printf("%d out of %d\n", test1, runs);
            System.out.println(test2 == runs);
            System.out.printf("%d out of %d\n", test2, runs);
            System.out.println(net.getNetwork());
        }
        {
            int runs = 10000;
            int test = 0;
            int test1 = 0;
            int test2 = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new DeleteReticulation(net);
                double logHR = op.propose();
                op.undo();
                test++;
                if(Networks.hasTheSameTopology(topo, net.getNetwork())) {
                    test1++;
                } else {
                    System.out.println(net.getNetwork().toString());
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
            int runs = 10000;
            int counter0 = 0, counter1 = 0, counter2 = 0, counter3 = 0, counter4 = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new AddReticulation(net);
                Network beforeProposal = Networks.readNetwork(net.getNetwork().toString());
                double logHR = op.propose();
                Networks.autoLabelNodes(net.getNetwork());
                if(!prior.isValid(net.getNetwork()) || !net.isValid() || !net.isUltrametric()) {
                    System.err.println("invalid add reticulation move!!!!");
                    System.exit(1);
                }
                while(net.getNetwork().getReticulationCount() != 3) {
                    op = new DeleteReticulation(net);
                    String prev = net.getNetwork().toString();
                    logHR = op.propose();
                    Networks.autoLabelNodes(net.getNetwork());
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
                        System.err.println(prev);
                        System.err.println(net.getNetwork().toString());
                        op.undo();
                    }
                    counter4++;
                    Networks.autoLabelNodes(net.getNetwork());

                    if(!prior.isValid(net.getNetwork()) || !net.isUltrametric() || !net.isValid()) {
                        System.err.println("!!!!!!!!! wrong Invaild move after undo \n" + beforeProposal.toString() + "\n"
                                + afterProposal + "\n"
                                + Networks.hasTheSameTopology(beforeProposal, net.getNetwork()));
                        System.err.println(net.getNetwork().toString());
                        System.exit(1);
                    }
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
            System.out.printf("%d out of %d valid network!!!\n", counter2, counter4);
        }
    }


}
