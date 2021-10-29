package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.test;
/*
 * @ClassName:   checkResultMCMC
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        10/7/21 9:38 PM
 */

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.*;


public class checkResultMCMC {
//    private Network topology = null;
//    private double proportion = 0.0;
    private Map<Network, Integer> network_list = new HashMap<>();
    private int burn_in;
    private int sample_frequency;
    private double lowest_ESS = 0;
    private double last_ESS = 0;

    /* Constructor */
//    public checkResultMCMC(String path) {
//        try {
//            BufferedReader br = new BufferedReader(new FileReader(path));
//            String s;
//
//            boolean begin = false;
//            boolean finished = false;
//            int cnt = 0;
//            while((s = br.readLine()) != null) {
//                s = s.trim();
//                if(s.contains("Logger:")) {
//                    begin = true;
//                }
//
//                if(begin && s.startsWith("[")) {
//                    cnt += 1;
//                    s = s.substring(s.indexOf("]") + 1);
//                    Network net = Networks.readNetwork(s);
//                    boolean inmap = false;
//                    for (Network n1: network_list.keySet()){
//                        if (Networks.hasTheSameTopology(n1, net)){
//                            network_list.put(n1, network_list.get(n1)+1);
//                            inmap = true;
//                            break;
//                        }
//                    }
//
//                    if (!inmap){
//                        network_list.put(net, 1);
//                    }
//                }
//
//
//            }
//
//            br.close();
//            Map<Network, Integer> res = sortByValue(network_list);
//            Map.Entry<Network, Integer> e = res.entrySet().iterator().next();
//            this.topology = e.getKey();
//            this.proportion = e.getValue()*1.0/cnt;
//
//        } catch(Exception e) {
////            e.printStackTrace();
//            System.err.println(path);
//        }
//
//    }


    public checkResultMCMC(int burnin, int sf){
        this.burn_in = burnin;
        this.sample_frequency = sf;
    }

    public checkResultMCMC(int burnin, int sf, double lowest_ESS){
        this(burnin, sf);
        this.lowest_ESS = lowest_ESS;
    }



    public int read_mcmc_out(String path){
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String s;
            int start = burn_in/sample_frequency + 1;
            double lastESS = 0;
            boolean begin = false;
            boolean finished = false;
            int cnt = 0;
            List<Network> net_list = new ArrayList<>();
            while((s = br.readLine()) != null) {
                s = s.trim();
                if(s.contains("Logger:")) {
                    begin = true;
                }
                if(begin && s.startsWith("[")) {
                    cnt += 1;
                    if (cnt >= start){
                        s = s.substring(s.indexOf("]") + 1);
                        Network net = Networks.readNetwork(s);
                        net_list.add(net);

                    }
                }
                else {
                    String ss[] = s.split("\\s+");
                    if(ss.length == 7 && Character.isDigit(ss[2].charAt(0)) && ss[2].charAt(ss[2].length() - 1) == ';'){
                        lastESS = Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));
                    }
                }
            }
            br.close();
            this.last_ESS = lastESS;
            if (lastESS > this.lowest_ESS){

                for (Network net :net_list){
                    boolean inmap = false;
                    for (Network n1: network_list.keySet()){
                        if (Networks.hasTheSameTopology(n1, net)){
                            network_list.put(n1, network_list.get(n1)+1);
                            inmap = true;
                            break;
                        }
                    }

                    if (!inmap){
                        network_list.put(net, 1);
                    }
                }
                return cnt - start + 1;
            }
            else{
                System.out.println(path+":not converged");
                return 0;
            }
        } catch(Exception e) {
            System.err.println(path);
            return -1;
        }
    }

    public Map<Network, Integer> get_network_list(){
        return sortByValue(network_list);
//        return this.network_list;
    }

    public void reset_network_list(){
        this.network_list = new HashMap<>();
    }

    public double getLast_ESS(){
        return last_ESS;
    }

    public StringBuilder summarize_multiple_chains(List<String> list_chains_path){
        int cnt = 0;
        StringBuilder sb = new StringBuilder();
        for ( String path : list_chains_path){
            int cnt_chain = read_mcmc_out(path);
            if (cnt_chain != -1){
                cnt += cnt_chain;
            }
        }
        sb.append("net,");
        sb.append("count,");
        sb.append("proportion\n");
        Map<Network, Integer> res = get_network_list();
//        Map<Network, Integer> res = sortByValue(network_list);
//        Map.Entry<Network, Integer> e = res.entrySet().iterator().next();
//        setTopology(e.getKey());
//        setProportion(e.getValue()*1.0/cnt);
        for (Network net: res.keySet()){
            System.out.println(net.toString()+"\t"+res.get(net)+"\t"+res.get(net)*1.0/cnt);
            sb.append("\"");
            sb.append(net.toString());
            sb.append("\",");
            sb.append(res.get(net));
            sb.append(",");
            sb.append(res.get(net)*1.0/cnt);
            sb.append("\n");


        }
        return sb;
    }

//    public Network getTopology(){
//        return topology;
//
//    }
//
//    public void setTopology(Network topology){
//        this.topology = topology;
//
//    }
//
//    public void setProportion(double proportion){
//        this.proportion = proportion;
//    }
//
//    public double getProportion(){
//        return proportion;
//    }



    public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
        List<Map.Entry<K, V>> list = new ArrayList<>(map.entrySet());
        list.sort(Map.Entry.comparingByValue());
        Collections.reverse(list);
        Map<K, V> result = new LinkedHashMap<>();
        for (Map.Entry<K, V> entry : list) {
            result.put(entry.getKey(), entry.getValue());
        }

        return result;
    }

    public static void summarize_chains_separate(){
        String directory = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/butterfly/seed3/";
        String output_path = directory;
        boolean murate = false;
        int burnin = 5000000;
        int sample_frequency = 5000;
        double lowest_ess = 200;
        StringBuilder sb = new StringBuilder();

        if (murate){
            output_path += "mcmc_murate_separate.csv";
        }
        else{
            output_path += "mcmc_original_separate.csv";
        }

        sb.append("seed,");
        sb.append("rank,");
        sb.append("net,");
        sb.append("count,");
        sb.append("proportion\n");

        for (int i = 11; i <= 20; i++){
            String chain = directory + i;
            if (murate){
                chain += "/murate/mcmc.out";
            }
            else{
                chain += "/original/mcmc.out";
            }

            checkResultMCMC cr = new checkResultMCMC(burnin, sample_frequency, lowest_ess);
            int cnt = cr.read_mcmc_out(chain);
            Map<Network, Integer> res = cr.get_network_list();

            int rank = 1;
            for (Network net : res.keySet()){
                sb.append(i);
                sb.append(",");
                sb.append(rank);
                sb.append(",");
                sb.append("\"");
                sb.append(net.toString());
                sb.append("\",");
                sb.append(res.get(net));
                sb.append(",");
                sb.append(res.get(net)*1.0/cnt);
                sb.append("\n");
                rank += 1;
            }

            cr.reset_network_list();

        }
        try (PrintWriter writer = new PrintWriter(new File(output_path))){
            writer.write(sb.toString());
        }catch (Exception e){
            System.err.println(e.getMessage());
        }
    }


    public static void summarize_chains_2_one(){
        String directory = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/butterfly/seed3/";
        String output_path = directory;
        boolean murate = true;
        int burnin = 5000000;
        int sample_frequency = 5000;
        int lowest_ess = 200;
        List<String> path_list = new ArrayList<>();
        if (murate){
            output_path += "mcmc_murate.csv";
        }
        else{
            output_path += "mcmc_original.csv";
        }
        for(int i = 11; i <= 20; i++){
            if (murate){
                path_list.add(directory+i+"/murate/mcmc.out");

            }
            else{
                path_list.add(directory+i+"/original/mcmc.out");

            }
        }
        checkResultMCMC cr = new checkResultMCMC(burnin, sample_frequency, lowest_ess);
        StringBuilder sb = cr.summarize_multiple_chains(path_list);
        try (PrintWriter writer = new PrintWriter(new File(output_path))){
            writer.write(sb.toString());
        }catch (Exception e){
            System.err.println(e.getMessage());
        }
    }


    public static void main(String[] args) {
//        summarize_chains_2_one();

        summarize_chains_separate();
    }

}
