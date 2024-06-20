package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.test;
/*
 * @ClassName:   test_rnBS
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        12/15/22 8:41 PM
 */

import com.google.common.collect.Lists;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.experiments.McmcseqPairwiseExperiments;
import edu.rice.cs.bioinfo.programs.phylonet.algos.summarize.majorTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.junit.Test;

import java.io.*;
import java.util.*;
import java.util.function.BiFunction;

public class MCMCrNBS {
    /* Constructor */


//    private String metric = "rNBS";
    private String metric = "Nakhleh";
//    private String metric = "NormRNBS";

    private Network _truenet = Networks.readNetwork("((((C:1.0)#H2:3.0::0.7)#H1:2.0::0.2,(D:3.0,#H2:2.0::0.3):3.0):2.0,(#H1:1.0::0.8,A:5.0):3.0);");
    List<String> _selectedLeaves;
    Map<String, String> _taxonMap = new HashMap<>();

    //    public static Map<String, BiFunction<double[][], double[][], Double>> metrics2 = new HashMap<>();
    public static Map<String, BiFunction<Network<BniNetwork>, Network<BniNetwork>, Double>> metrics = new HashMap<>();

    static {
//        metrics2.put("WAPD", McmcseqPairwiseExperiments::getWAPD);
//        metrics2.put("normWAPD", McmcseqPairwiseExperiments::getNormWAPD);
//        metrics2.put("APD", McmcseqPairwiseExperiments::getAPD);
//        metrics2.put("normAPD", McmcseqPairwiseExperiments::getNormAPD);
        metrics.put("rNBS", McmcseqPairwiseExperiments::getRNBS);
        metrics.put("NormRNBS", McmcseqPairwiseExperiments::getNormRNBS);
        metrics.put("Nakhleh", McmcseqPairwiseExperiments::getNakhleh);
    }

    public MCMCrNBS() {
        initialize();

    }
    public MCMCrNBS(List<String> selectedLeaves) {
        _selectedLeaves = selectedLeaves;
        initialize();

    }


    public MCMCrNBS(List<String> selectedLeaves, Map<String, String> taxonmap) {
        _selectedLeaves = selectedLeaves;
        _taxonMap = taxonmap;
        initialize();

    }

    private void initialize(){
        String truenetstr = "(((A:5,((B:2,(C:1)#H1:1::0.7):2)#H2:1::0.8):3,((#H1:2::0.3,D:3):3,#H2:2::0.2):2):2,E:10):0;";
//        String truenetstr = "(((A_0:5,((B_0:2,(C_0:1)#H1:1::0.7):2)#H2:1::0.8):3,((#H1:2::0.3,D_0:3):3,#H2:2::0.2):2):2,E_0:10):0;";
        Network truenet = Networks.readNetwork(truenetstr);
        double theta = 0.04;

        if (_selectedLeaves == null || _selectedLeaves.isEmpty()){
            _truenet = truenet;
        }
        else{
            Tuple<Network, Map<NetNode, NetNode>> tuple = SuperNetwork3.getSubNetwork(truenet, _selectedLeaves, true);
            _truenet = tuple.Item1;
        }

        if (!_taxonMap.isEmpty()){
            for(Object o: _truenet.getLeaves()){
                NetNode n = (NetNode) o;
                n.setName(_taxonMap.get(n.getName()));
            }
        }
        Networks.scaleNetwork(_truenet, theta/2);
        System.out.println(_truenet.toString());
    }


    public void calcDistMCMC(String filename, int start){
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String s;
            List<Double> dist = new ArrayList<>();
            boolean begin = false;
            int i = 0;
            while((s = in.readLine()) != null) {

                if(begin) {
                    i++;
                    if (i >= start){
                        if (s.startsWith("[")) {
                            Network curSample = Networks.readNetworkWithRootPop(s);

//                            System.out.println(curSample.getRoot().getRootPopSize());
//                            double theta = curSample.getRoot().getRootPopSize();
//                            System.out.println(s);
//                            System.out.println(2/theta);
//                            Networks.scaleNetwork(curSample, 2/theta);
//                            System.out.println(curSample.toString());
//                            System.out.println(_truenet.toString());
                            dist.add(metrics.get(metric).apply(_truenet, curSample));

                        }
                    }

                } else {
                    if(s.contains("Logger")) {
                        begin = true;
                    }
                }
            }

            in.close();
            System.out.println(dist);

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public double computeDist(Network net, String metric){

        System.out.println(_truenet.toString());
        System.out.println(net.toString());
        return metrics.get(metric).apply(_truenet, net);
    }

    public static void calcEnum(){
        String path = "";
        Map<String, List> file2samples = new HashMap<>();
        List<Double> sampled_dist = new ArrayList<>();
        List<String> selectedLeaves = Arrays.asList(new String[] {"A", "C", "D"});

        for(int i = 1; i <= 288; i++) {
            String filename = path+"/.out";
            MCMCrNBS bs = new MCMCrNBS(selectedLeaves);
            bs.calcDistMCMC(filename, 1);

        }
    }

    public static void calcEnumres(){
        String num = "9";
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/simulation/taxa5_1individual/enumerate/enum/mcmc"+num+"/res.csv";
        String output_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/simulation/taxa5_1individual/enumerate/enum/mcmc"+num+"/res_dist.csv";
//        List<String> selectedLeaves = Arrays.asList(new String[] {"A", "D", "E"});
//        String tm = "1:A; 2:D; 3:E";

//        List<String> selectedLeaves = Arrays.asList(new String[] {"A", "C", "D"});
//        String tm = "1:A; 2:D; 3:C";

        List<String> selectedLeaves = Arrays.asList(new String[] {"B", "D", "E"});
        String tm = "1:B; 2:D; 3:E";
        String [] tmarr = tm.split(";");
        Map<String, String> taxaMap = new HashMap<>();
        for (String pair: tmarr){
            String [] ab = pair.trim().split(":");
            taxaMap.put(ab[1], ab[0]);
        }


        MCMCrNBS bs = new MCMCrNBS(selectedLeaves, taxaMap);
        String[] metrics = {"Nakhleh", "rNBS", "NormRNBS"};
        StringBuilder output = new StringBuilder();
        try {
            BufferedReader in = new BufferedReader(new FileReader(path));
            String s;
            List<Double> distList = new ArrayList<>();
            boolean begin = false;
            int i = 0;
            while ((s = in.readLine()) != null) {
//                System.out.println(s);
                if (begin){
                    String[] arr = s.split("\"");
                    Network net = Networks.readNetwork(arr[1]);
                    output.append(s);
                    output.append(",");
                    for(String metric : metrics){
                        double dist = bs.computeDist(net, metric);
                    System.out.println(dist);
                        output.append(dist);
                        output.append(",");
                    }
                    output.append("\n");

                }
                else{
                    output.append(s);
                    output.append(",");
                    for(String metric : metrics){
                        output.append(metric);
                        output.append(",");
                    }
                    output.append("\n");
                }
                begin = true;

            }
            System.out.println(output);
            BufferedWriter outAll = new BufferedWriter(new FileWriter(output_path));

//            PrintWriter out = new PrintWriter(output_path);
            outAll.flush();
            outAll.write(output.toString());
            outAll.close();

        }catch (Exception e){
            System.err.println(e);
        }

    }

    public static void MCMCDistance(){
//        List<String> selectedLeaves = Arrays.asList(new String[] {"A", "D", "E"});
//        List<String> selectedLeaves = Arrays.asList(new String[] {"A", "C", "D"});
//        List<String> selectedLeaves = Arrays.asList(new String[] {"A", "B", "C", "D", "E"});

//        List<String> selectedLeaves = Arrays.asList(new String[] {"B", "D", "E"});

//
//        String filename = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/simulation/taxa5_1individual/heter/mcmc6/mcmc_sgt.out";
        String filename = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/simulation/taxa5_1individual/original_op/mcmc/mcmc_sgt.out";
//        String filename = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/simulation/taxa5_1individual/heter/mcmc9/mcmc_sgt.out";

//        MCMCrNBS bs = new MCMCrNBS(selectedLeaves);
        MCMCrNBS bs = new MCMCrNBS();

        bs.calcDistMCMC(filename, 1);
    }

    public static void main(String[] args) {
//        MCMCDistance();

//        double x = metrics.get("rNBS").apply(Networks.readNetwork("((D:10.3398810817137994,C:2.3398810817137994):0.17816904346420231,A:0.5180501251780016);"), Networks.readNetwork("((D:0.3398810817137994,C:0.3398810817137994):0.17816904346420231,A:0.5180501251780016);"));
//        System.out.println(x);
//        calcEnumres();
        Network n1 = Networks.readNetwork("(3:0.1717052634863686,((1:0.06400114979146808)I2#H1:0.08021684469799646::0.45689811700722466,(I2#H1:0.06388604553294634::0.5431018829927754,2:0.12788719532441442)I3:0.01633079916505012)I1:0.02748726899690407)I0;\n");
        Network n2 = Networks.readNetwork("((((2:0.07046803098170829)I5#H2:0.06999617872142991::0.6443831228683109,1:0.1404642097031382)I3:0.01510798163628957)I1#H1:0.07561150253708232::0.22621238703772623,(3:0.18229306450171334,(I5#H2:0.0974957148609727::0.3556168771316891,I1#H1:0.012391554503253216::0.7737876129622738)I2:0.014329318659032364)I4:0.04889062937479674)I0;\n");
        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(n1, n2);
        System.out.println(closest.Item1.getReticulationCount());
        System.out.println(closest.Item2.getReticulationCount());
        List<NetNode> netnodeList = Lists.newArrayList(n2.getNetworkNodes());

        majorTree.removeReti(n2, netnodeList, 0);
        System.out.println(n2.toString());
        System.out.println(closest.Item3);
    }
}
