package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.param;

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
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 3/4/16.
 * Scales the time of the root by a random factor.
 */
public class ScaleRootTime extends NetworkOperator {

    private double _scaleFactor = 0.90;
    private double _upperLimit = 1.0 - 1e-6;
    private double _lowerLimit = 1e-6;

    public double scale;
    private double _logHR;
    private Double _oldHeight;

    public ScaleRootTime(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        // reset
        _oldHeight = null;

        scale = getScaler(); // scale > 1 => increase height => may violate

        NetNode<NetNodeInfo> root = _network.getNetwork().getRoot();
        _oldHeight = root.getData().getHeight();
        double newHeight = _oldHeight * scale;
        double[] bounds = _network.getLowerAndUpperBound(root);

        if(newHeight < bounds[0] || newHeight > Utils.NET_MAX_HEIGHT) {
            _logHR = Utils.INVALID_MOVE;
            _violate = false;
        } else {
            setNodeHeight(root, newHeight);
            _logHR = -Math.log(scale);
            _violate = scale > 1.0;
        }
        return _logHR;
    }

    @Override
    public void undo() {
        if (_logHR == Utils.INVALID_MOVE) return;
        setNodeHeight(_network.getNetwork().getRoot(), _oldHeight);
    }

    @Override
    public String getName() {
        return "Scale-Root-Time";
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
        {
            UltrametricNetwork net = new UltrametricNetwork("((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0:3.0;"
                                                           , geneTrees);
            System.out.println(geneTrees.get(0).getTree().toNewick());
            System.out.println(net.getNetwork().toString());
            int runs = 10000;
            int test = 0;
            int test2 = 0;
            for(int i = 0; i < runs; i++) {
                ScaleRootTime op = new ScaleRootTime(net);
                double logHR = op.propose();
                if(net.isValid()) {
                    test2++;
                } else {
                    op.undo();
                }
                if(net.isUltrametric()) test++;
            }
            System.out.println(net.getNetwork());
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
            System.out.printf("%d out of %d\n", test2, runs);
        }
        {
            UltrametricNetwork net = new UltrametricNetwork("((((A:0.5,(C:0.5,(G:0.424733668924115)I10#H3:0.27526633107588505::0.95)I6:1.2148953368828783)I5:2.0902392667517367,(((R:0.23468712)I8#H2:0.57719209683010764::0.76)I2#H1:0.4723366810412248::0.29,(Q:0.5,I10#H3:0.24368001872495645::0.05)I4:1.4561500003545642)I9:2.081453711192161)I3:0.8324858688602822,(L:0.30183286,I8#H2:0.8699237693946538::0.24)I7:0.7369103845180647)I1:1.427856628628852,I2#H1:0.8195085232672825::0.71)I0:3.0;"
                    , geneTrees);
            System.out.println(net.getNetwork());
            List<Double> heights = net.getOldHeights();
            int runs = 10000;
            int test = 0;
            int test2 = 0;
            for(int i = 0; i < runs; i++) {
                NetworkOperator op = new ScaleRootTime(net);
                double logHR = op.propose();
                op.undo();
                int idx = 0;
                for(NetNode<NetNodeInfo> node : Networks.postTraversal(net.getNetwork())) {
                    if(Math.abs(node.getData().getHeight() - heights.get(idx++)) > 0.000001) {
                        System.out.println(node.getData().getHeight() + " vs " + heights.get(idx-1));
                        test--;
                        break;
                    }
                }
                test++;
                if(net.isUltrametric()) test2++;
            }
            System.out.println(test == runs);
            System.out.printf("%d out of %d\n", test, runs);
            System.out.println(test2 == runs);
            System.out.printf("%d out of %d\n", test2, runs);
            System.out.println(net.getNetwork());
        }
    }

}
