package edu.rice.cs.bioinfo.programs.phylonet.algos.butterfly;
/*
 * @ClassName:   TrinetStats
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        2/18/22 8:21 PM
 */


import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.util.*;

public class TrinetStats {

    public static boolean printDetails_ = false;


    private static class ValueComparator implements Comparator<Map.Entry<String, Tuple<Double, Double>>> {
        public int compare(Map.Entry<String, Tuple<Double, Double>> mp1, Map.Entry<String, Tuple<Double, Double>> mp2)
        {
            if (mp2.getValue().Item2 - mp1.getValue().Item2 > 0.00001){
                return 1;
            }
            else if (mp2.getValue().Item2 - mp1.getValue().Item2 < -0.00001){
                return -1;

            }
            else{
                return 0;
            }
        }
    }

    public static StringBuilder analyze1Network(String filename, int start, Map<String, List<String>> file2samples, List<Double> top_frequency_list, boolean report_ML){
        try {
            double priorESS = 0.0;
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
            StringBuilder sb = new StringBuilder();
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
                        if(ss.length == 9 && Character.isDigit(ss[2].charAt(0)) && ss[2].charAt(ss[2].length() - 1) == ';'){
                            lastESS = Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));
                            priorESS = Double.parseDouble(ss[5].substring(0, ss[5].length() - 1));
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
                return null;
            }

            Network<Object> topNet = Networks.readNetworkWithRootPop(topSample);
            List<String> selectedLeaves = new ArrayList<>();

            for(NetNode leaf : topNet.getLeaves()) {
                selectedLeaves.add(leaf.getName());
            }
            top_frequency_list.add(maxValue*1.0/totalValue);

//            System.out.print("\""+Networks.readNetwork(MLSample).toString() + "\",");
//            System.out.print(Networks.readNetwork(MLSample).getReticulationCount()+",");
//            System.out.print("\""+topNet.toString()+"\",");
//            System.out.print(topNet.getReticulationCount()+",");
//            System.out.print(lastESS);
//            System.out.print("\n");
            sb.append("\""+Networks.readNetwork(MLSample).toString() + "\",");
            sb.append(Networks.readNetwork(MLSample).getReticulationCount()+",");
            sb.append("\""+topNet.toString()+"\",");
            sb.append(topNet.getReticulationCount()+",");
            sb.append(lastESS);
            sb.append("\n");
            return sb;
//            if(report_ML){
//                System.out.print(Networks.readNetwork(MLSample).toString() + ",");
//                System.out.print(Networks.readNetwork(MLSample).getReticulationCount()+"\n");
////                if (Networks.readNetwork(MLSample).getReticulationCount() >= 1){
//////                    System.out.println(filename);
//////                    System.out.println("reti in ML ="+Networks.readNetwork(MLSample).getReticulationCount()+": "+selectedLeaves);
////                    System.out.println(Networks.readNetwork(MLSample).toString());
////                }
//            }
//            else{
//                if (topNet.getReticulationCount() >= 1){
//                    System.out.println(filename);
//                    System.out.println("reti in ML ="+topNet.getReticulationCount()+": "+selectedLeaves);
//                    System.out.println(topNet.toString());
//                }
//            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    public static void analyzeNetwork(String resultFolder) {
//        int chainlen = 60000000;
        int chainlen = 30000000;
        int burnin = 1000000;
        int sample_freq = 5000;
        List<List<String>> topsample_reti = new ArrayList<>();
        List<List<String>> MLsample_reti = new ArrayList<>();
        List<Double> top_frequency_list = new ArrayList<>();
        int start = (int)(burnin / sample_freq) + 1;
        int end = (int)(chainlen / sample_freq);
        double priorESS = 0.0;
        File path = new File(resultFolder);
        boolean report_ML = false;

        List<String> filenames = new ArrayList<>();

        File[] files = path.listFiles();
//        for (int i = 0; i < files.length; i++){
//            if (files[i].isFile()){ //this line weeds out other directories/folders
//                if(files[i].toString().endsWith(".out")) {
//                    //System.out.println(files[i]);
//                    filenames.add(files[i].toString());
//                }
//            }
//        }
//        for (int i = 301; i <= 560; i++) {
        for(int i = 5556251; i <=5556307; i++){
//            String filename = resultFolder + "/slurm-" + String.valueOf(i) + ".out";
            String filename = resultFolder + "/slurm-" + String.valueOf(i) + "_0.out";
            if  (new File(filename).exists()){
                filenames.add(filename);
            }

            filename = resultFolder + "/slurm-" + String.valueOf(i) + "_1.out";
            if  (new File(filename).exists()){
                filenames.add(filename);
            }
        }



//        Collections.sort(filenames);

        Map<String, List<String>> file2samples = new HashMap<>();

        try{
            BufferedWriter writer = new BufferedWriter(new FileWriter(resultFolder+"/res.csv", true));
//            writer.append(' ');
//            writer.append(str);


            for(String filename : filenames) {
                file2samples.put(filename, new ArrayList<>());
                System.gc();
                System.out.print(filename+",");
                writer.write(filename+",");
                StringBuilder sb = analyze1Network(filename, start, file2samples, top_frequency_list, report_ML);
                writer.write(sb.toString());
            }
            writer.close();
        }catch (IOException e){
            System.err.println(e);
        }

//        for (double x :top_frequency_list){
//            System.out.println(x);
//        }
    }

    public static void trinets(){
        String filename = "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/iqtree/sig_triplets_exact_bonferroni_map.csv";
        String slurm_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/trinet2merge_500/slurm-";
        int chainlen = 80000000;
        int burnin = 1000000;
        int sample_freq = 5000;
        List<Double> top_frequency_list = new ArrayList<>();
        boolean report_ML = true;
        try{
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String s = null;
            int start = (int)(burnin / sample_freq) + 1;
            while((s = in.readLine()) != null) {
                String[] arr = s.split(",");
                String slurm_path = slurm_dir + arr[0] + ".out";
                Map<String, List<String>> file2samples = new HashMap<>();
                file2samples.put(slurm_path, new ArrayList<>());
                analyze1Network(slurm_path, start, file2samples, top_frequency_list, report_ML);
            }
        }
        catch (Exception e){
            System.out.println("exception");
        }
    }


    public static void main(String[] args) {
//        String folder = "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/trinet2merge/";
//        String folder = "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/trinet2merge_500/";
        String folder = "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/trinets13/invgamma/";
        analyzeNetwork(folder);
//        trinets();
    }
}
