package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.InferTreeWrapper;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 3/23/18
 * Time: 5:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimSeqInNetworkByMSSeqGen {
    static public List<Map<String, String>> generateSeqences(Network network, Map<String, List<String>> species2alleles, int numGTs, int lenPerGT, double theta, double freq[], double rates[], String MSPath, String SeqGenPath, String gtFilePath) {
        SimGTInNetworkByMS simGT = new SimGTInNetworkByMS();
        List<Tree> gts = simGT.generateGTs(network, species2alleles, numGTs, MSPath);
        List<Map<String, String>> lociSeq = new ArrayList<>();

        List<List<MutableTuple<Tree, Double>>> treeList = new ArrayList<>();
        for(Tree tree : gts) {
            /*SuperNetwork superNetwork = new SuperNetwork(new ArrayList<>());
            superNetwork.setTrueNetwork(Networks.readNetwork(tree.toNewick()));
            List<String> alleles = new ArrayList<>();
            alleles.add("A_0");
            alleles.add("A_1");
            alleles.add("B_0");
            alleles.add("B_1");
            alleles.add("D_0");
            alleles.add("D_1");

            Network net = superNetwork.getSubNetwork(Networks.readNetwork(tree.toNewick()), alleles, true, false);
            tree = Trees.readTree(net.toString());*/
            treeList.add(Arrays.asList(new MutableTuple<Tree, Double>(Trees.readTree(tree.toNewick()), 1.0)));

        }

        InferTreeWrapper inferTreeWrapper = new InferTreeWrapper();
        List<MutableTuple<Tree, Double>> distinctTrees = new ArrayList<>();
        inferTreeWrapper.summarizeData(treeList, null, distinctTrees);


        for(MutableTuple<Tree, Double> gtTuple : distinctTrees) {
            System.out.println(gtTuple.Item1.toNewick() + " " + gtTuple.Item2);
        }

        for(Tree gt : gts) {
            /*SuperNetwork superNetwork = new SuperNetwork(new ArrayList<>());
            superNetwork.setTrueNetwork(Networks.readNetwork(gt.toNewick()));
            List<String> alleles = new ArrayList<>();
            alleles.add("A_0");
            alleles.add("A_1");
            alleles.add("B_0");
            alleles.add("B_1");
            alleles.add("D_0");
            alleles.add("D_1");
            Network net = superNetwork.getSubNetwork(Networks.readNetwork(gt.toNewick()), alleles, true, false);
            gt = Trees.readTree(net.toString());
            */

            lociSeq.add(SimSeqInGTBySeqGen.execute(gt, theta, freq, rates, lenPerGT, null, SeqGenPath));
        }

        if(gtFilePath != null) {
            try {
                PrintWriter out = new PrintWriter(gtFilePath);
                for(Tree gt : gts) {
                    out.println(gt.toNewick());
                }
                out.close();
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
        }

        return lociSeq;
    }

    public static void main(String []args) {
        Network network = Networks.readNetwork("(((((C:0.2)#H1:0.6::0.5,B:0.8):0.4)#H2:0.2::0.5,A:1.4):0.6,((D:0.4,#H1:0.2::0.5):1.2,#H2:0.4::0.5):0.4);");
        String mspath = "/Users/zhujiafan/Documents/Luay/msdir/ms";
        String seqgenpath = "/Users/zhujiafan/Documents/Luay/Seq-Gen.v1.3.3/source/seq-gen";
        double[] base = {0.2112, 0.2888, 0.2896, 0.2104};
        double[] trans = {0.2173, 0.9798, 0.2575, 0.1038, 1, 0.2070};
        List<Map<String, String>> lociSeq = SimSeqInNetworkByMSSeqGen.generateSeqences(network, null, 10, 50, 0.01, base, trans, mspath, seqgenpath, null);
        System.out.println("Done");
    }
}
