package edu.rice.cs.bioinfo.programs.phylonet.algos.butterfly;
/*
 * @ClassName:   TestMerger
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        2/21/22 2:41 PM
 */

import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class TestMerger {
    /* Constructor */
    public TestMerger() {
    }

    public static void testPipelineStage2_1(){
        String resultFolder = "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/inferred_trinets";
        int chainlen = 20000000;
        int burnin = 2000000;
        int sample_freq = 5000;
        SNOptions options = new SNOptions();
        options.outgroup = "Avan";
        options.trustReticulationTime = false;
        options.reconcileHeights = true;

        File path = new File(resultFolder);
        List<String> filenames = new ArrayList<>();

        File[] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].toString().endsWith(".out")) {
                    //System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);
        SNSummary snSummary = Pipeline.stage2_1(filenames, chainlen, burnin, sample_freq, options);
        Network inferred = snSummary.inferredNetwork;
        System.out.println("Final network:");
        System.out.println(inferred.toString());
        System.out.println("end");
    }

    public static Network trinetMerge(List<Network> trinetlist, SNOptions options){
        System.out.println();
        SNProblem problem = new SNProblem();
        Random random = new Random(12345678L);
        for(Network network : trinetlist) {
            Network newnet = network.clone();
            NetworkUtils.alterHeights(newnet, random);
            problem.AddSubNetwork(newnet, "", 1.0);
        }

        SNSummary summary = SNSolver.Solve(problem, options);
        return summary.inferredNetwork;
    }

    public static List<Network> getTrinets(String resultFolder, int trinetType) {
        int chainlen = 20000000;
        int burnin = 1000000;
        int sample_freq = 5000;
        List<Network> trinetList = new ArrayList<>();
        int start = (int)(burnin / sample_freq) + 1;
        int end = (int)(chainlen / sample_freq);

        File path = new File(resultFolder);

        List<String> filenames = new ArrayList<>();

        File[] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].toString().endsWith(".out")) {
                    //System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);

        Map<String, List<String>> file2samples = new HashMap<>();

        for(String filename : filenames) {
            file2samples.put(filename, new ArrayList<>());

            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                int index = 0;
                String topSample = null;
                String MLSample =  null;
                double maxML = Double.NEGATIVE_INFINITY;
                double MLvalue = 0;
                Map<Network, Integer> count = new HashMap<>();
                boolean begin = false;
                int total = 0;
                double lastESS = 0.0;
                boolean updateMLsample = false;

                while((s = in.readLine()) != null) {

                    if(begin) {
                        if (s.startsWith("[")) {
                            Network curSample = Networks.readNetworkWithRootPop(s);
                            if (updateMLsample){
                                MLSample = s;
                            }
                            total++;
                            if(total >= start) {
                                file2samples.get(filename).add(s);
                                boolean exist = false;
                                for (Network net : count.keySet()) {
                                    if (Networks.hasTheSameTopology(net, curSample)) {
                                        count.put(net, count.get(net) + 1);
                                        exist = true;
                                        break;
                                    }
                                }
                                if (!exist) {
                                    count.put(curSample, 1);
                                }
                            }
                        } else {
                            String ss[] = s.split("\\s+");
                            if(ss.length == 7 && Character.isDigit(ss[2].charAt(0)) && ss[2].charAt(ss[2].length() - 1) == ';'){
                                lastESS = Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));
                                MLvalue = Double.parseDouble(ss[1].substring(0, ss[1].length() - 1));
                                if (MLvalue > maxML){
                                    maxML = MLvalue;
                                    updateMLsample = true;
                                }
                                else{
                                    updateMLsample = false;
                                }
                            }


                        }
                    } else {
                        if(s.contains("Logger")) {
                            begin = true;
                        }
                    }
                    index++;
                }

                int totalValue = 0;
                int maxValue = 0;
                for(Network net : count.keySet()) {
                    if(maxValue < count.get(net)) {
                        topSample = Networks.getFullString(net);
                        maxValue = count.get(net);
                    }
                    totalValue += count.get(net);
                }
                in.close();

                if(topSample == null) {
                    continue;
                }

                Network<Object> topNet = Networks.readNetworkWithRootPop(topSample);

                if (trinetType == 1){
                    trinetList.add(Networks.readNetwork(topSample));
                }
                else if(trinetType == 2){
                    trinetList.add(Networks.readNetwork(MLSample));
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return trinetList;
    }


    public static void merge(){
        String resultFolder = "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/inferred_trinets";
        List<Network> trinetlist = getTrinets(resultFolder, 2);
        SNOptions options = new SNOptions();
        options.outgroup = "Avan";
        options.reconcileHeights = true;
        options.trustReticulationTime = false;

        Network net = trinetMerge(trinetlist, options);
        System.out.println("Inferred net:");
        System.out.println(net);
        System.out.println(net.toString());

    }


    public static void main(String[] args) {
        testPipelineStage2_1();
//        merge();
    }
}
