package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.test;
/*
 * @ClassName:   Simulations
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        12/20/22 3:47 PM
 */

import com.google.common.collect.Multimap;
import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.RootedNetworkBranchScore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.experiments.McmcseqPairwiseExperiments;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.util.CombinatoricsUtils;

import java.io.*;
import java.util.*;

public class Simulations {

    List<String> sampleNames = new ArrayList<>();
    Multimap<String, String> species2samples = null;
    Map<String, List<String>> species2alleles = new HashMap<>();
    Map<String, String> alleles2species = new HashMap<>();
    Map<String, Map<String, String>> loci2seq = new HashMap<>();
    int nchar = 0;
    boolean sgt = true;
    String outgroup = "Lampropholis_guichenoti";
    String gtoutgroup = "SP07_indexing28_h0";
    String root;

    /* Constructor */
    public Simulations() {

    }

    static List<String> readIQTree(String iqpath) {

//        List<String> filenames = new ArrayList<>();
//        Map<String, String> filename2locusname = new HashMap<>();
//        Map<Integer, String> locus2treestring = new TreeMap<>();
        List<String> treelist = new ArrayList<>();


        try {
            BufferedReader in = new BufferedReader(new FileReader(iqpath));
            String s;
            int i = 1;
            while((s = in.readLine()) != null){
                Tree tree = Trees.readTree(s);
                treelist.add( s);
                i++;
            }

        } catch (IOException e) {
            e.printStackTrace();
        }


        return treelist;
    }


    public static void checkSortedEnumNet(){
        Network truenet = Networks.readNetwork("((((3:1.0)#H2:3.0::0.7)#H1:2.0::0.2,(2:3.0,#H2:2.0::0.3):3.0):2.0,(#H1:1.0::0.8,1:5.0):3.0);");
        Networks.scaleNetwork(truenet, 0.01);
//        Network inferrednet0 = Networks.readNetwork("(((3:0.022579202994172872)I2#H1:0.02401082229252763::0.40929675535218596,2:0.046590025286700504)I1:0.06605400980489412,(I2#H1:0.048184331500900596::0.590703244647814,1:0.07076353449507347)I3:0.04188050059652115)I0;");
//        Network inferrednet1 = Networks.readNetwork("(((3:0.01564171225815039)I1#H2:0.01576580185875656::0.4432980705462146)I2#H1:2.1136118496985343::0.0697042997040119,((I2#H1:0.015816313878307946::0.9302957002959881,2:0.04722382799521489)I4:0.06325513778809005,(1:0.06632485389956679,I1#H2:0.05068314164141639::0.5567019294537854)I5:0.04415411188373816)I3:2.0345403980321364)I0;");

//        System.out.println(Networks.computeDistanceBetweenTwoNetworks(truenet, inferrednet0));
//        System.out.println(Networks.computeDistanceBetweenTwoNetworks(truenet, inferrednet1));

        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/simulation/taxa5_1individual/heter/mcmc4/sorted_post.txt";
//        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/simulation/taxa5_1individual/heter/mcmc4/sorted.txt";

        List<Network> networklist = new ArrayList<>();
        try {
            BufferedReader in = new BufferedReader(new FileReader(path));
            String s;
            int i = 1;
            while((s = in.readLine()) != null){
                Network net = Networks.readNetwork(s);
                double dist = Networks.computeDistanceBetweenTwoNetworks(net, truenet);
                double dist_rnbs = RootedNetworkBranchScore.compute(net, truenet);
//                System.out.println(dist+"-"+net);
                System.out.println(dist_rnbs +"-" + net);
                networklist.add(net);
                i++;
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void checkEnumNet() {
        Network inferred = Networks.readNetwork("(3:0.17631349459062517,(1:0.14498699154422978,2:0.14498699154422978)I1:0.03132650304639539)I0;");
        Network truenet = Networks.readNetwork("(3:10.0,(2:8.0,1:8.0):2.0);");
        Networks.scaleNetwork(truenet, 0.01);
        double dist_rnbs = RootedNetworkBranchScore.compute(inferred, truenet);
        System.out.println(dist_rnbs);
    }

    public static void main(String[] args) {
//        checkSortedEnumNet();



    }
}
