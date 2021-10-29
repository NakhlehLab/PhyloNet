package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.test;
/*
 * @ClassName:   SummarizeReplica
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        10/11/21 11:30 AM
 */

import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SummarizeReplica {
//    Network truenet;
    StringBuilder sb = new StringBuilder();
    Map<String, Network> truenetmap = new HashMap<>();

    /* Constructor */
    public SummarizeReplica() {
        truenetmap.put("tree", Networks.readNetwork("(G:5.0,(L:2.0,(Q:1.0,A:1.0):1.0):3.0);"));
        truenetmap.put("net", Networks.readNetwork("(((Q:0.5)I8#H1:3.5::0.3,(L:2.0,(I8#H1:0.5::0.7,A:1.0)I4:1.0)I3:2.0)I2:1.0,G:5.0)I0;"));
        truenetmap.put("AZ", Networks.readNetwork("((((A:0.1,B:0.1):0.2,C:0.3):0.1,D:0.4):10,E:10.4);"));
    }

    public void set_title_mcmc(){
        sb.append("type,");
        sb.append("index,");
        sb.append("murate,");
        sb.append("homo,");
        sb.append("num,");
        sb.append("net,");
        sb.append("proportion,");
        sb.append("ESS,");
        sb.append("distance,");
        sb.append("closet_distance,");
        sb.append("inferred_reti,");
        sb.append("true_reti,");
        sb.append("inferred_backbone_reti,");
        sb.append("true_backbone_reti,");

        sb.append("\n");
    }


    public void summarize_mcmc(String dir, String type, boolean murate, boolean homo, int burnin, int sf){
        Network true_net = truenetmap.get(type);
        int correct_cnt = 0;

        for (int i = 1; i <= 10; i++){
            String path = dir+"/"+type+"/"+String.valueOf(i)+"/";
            if (homo){
                path += "homo/mcmc/";
            }
            else{
                path += "heter_locus/mcmc/";
            }
            if (murate){
                path += "murate/";
            }
            path += "mcmc.out";

            checkResultMCMC cr = new checkResultMCMC(burnin, sf);
            int cnt = cr.read_mcmc_out(path);
            System.out.println(path);
            Map<Network,Integer> net_list = cr.get_network_list();
            Map.Entry<Network, Integer> e = net_list.entrySet().iterator().next();
            Network net = e.getKey();
            sb.append(type);
            sb.append(",");
            sb.append(i);
            sb.append(",");
            sb.append(murate);
            sb.append(",");
            sb.append(homo);
            sb.append(",");
            sb.append(i);
            sb.append(",");
            sb.append("\"");

            if(net != null){
                sb.append(net.toString());
                sb.append("\"");
                sb.append(",");
                sb.append(e.getValue()*1.0/cnt);
                sb.append(",");
                sb.append(cr.getLast_ESS());
                sb.append(",");
                sb.append(Networks.computeDistanceBetweenTwoNetworks(net, true_net));
                sb.append(",");


                Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(net, true_net);
                sb.append(closest.Item3);
                sb.append(",");

                sb.append(net.getReticulationCount());
                sb.append(",");

                sb.append(true_net.getReticulationCount());
                sb.append(",");

                sb.append(closest.Item1.getReticulationCount());
                sb.append(",");

                sb.append(closest.Item2.getReticulationCount());
            }
            else{
                System.err.println(path);
                sb.append("\"");
                sb.append(",");
                sb.append(",");
                sb.append(",");
                sb.append(",");

                sb.append(",");

                sb.append(",");

                sb.append(",");

            }

            sb.append("\n");

        }

    }


    public void set_title_ML(){
        sb.append("type,");
        sb.append("index,");
        sb.append("homo,");
        sb.append("gtt,");
        sb.append("reti_set,");
        sb.append("net,");
        sb.append("likelihood,");

        sb.append("distance,");
        sb.append("closet_distance,");
        sb.append("inferred_reti,");
        sb.append("true_reti,");
        sb.append("inferred_backbone_reti,");
        sb.append("true_backbone_reti,");

        sb.append("\n");
    }

    public void summarize_ML(String dir, String type, String gtt, boolean homo){
        Network true_net = truenetmap.get(type);

        for (int i = 1; i <= 10; i++){
            String path = dir+"/"+type+"/"+String.valueOf(i)+"/";
            if (homo){
                path += "homo/ML/";
            }
            else{
                path += "heter_locus/ML/";
            }
            path += gtt+"/ML_";
            for (int reti=0; reti <=3; reti++){
                String mlpath = path + reti + ".out";

                checkResultML cr = new checkResultML(mlpath);
                Network net = cr.getTopology();

                sb.append(type);
                sb.append(",");
                sb.append(i);
                sb.append(",");
                sb.append(homo);
                sb.append(",");
                sb.append(gtt.equals("truegt"));
                sb.append(",");
                sb.append(reti);
                sb.append(",");
                sb.append("\"");
                if (net != null){
                    sb.append(net.toString());
                }
                else{
                    System.out.println(mlpath);
//                    sb.append(" ");
                }

                sb.append("\"");
                sb.append(",");
                sb.append(cr.getLikelihood());
                sb.append(",");
                if (net != null){
                    sb.append(Networks.computeDistanceBetweenTwoNetworks(net, true_net));
                    sb.append(",");


                    Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(net, true_net);
                    sb.append(closest.Item3);
                    sb.append(",");

                    sb.append(net.getReticulationCount());
                    sb.append(",");

                    sb.append(true_net.getReticulationCount());
                    sb.append(",");

                    sb.append(closest.Item1.getReticulationCount());
                    sb.append(",");

                    sb.append(closest.Item2.getReticulationCount());

                }
                else{
//                    sb.append(" ");
                    sb.append(",");
//                    sb.append(" ");
                    sb.append(",");
//                    sb.append(" ");
                    sb.append(",");
//                    sb.append(" ");
                    sb.append(",");
//                    sb.append(" ");
                    sb.append(",");
//                    sb.append(" ");
                }

                sb.append("\n");

            }
        }
    }


    public static void checkTopologies(String output_path, String method, int burnin, int sf){
        try (PrintWriter writer = new PrintWriter(new File(output_path))) {
            String dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/replica2/";
            SummarizeReplica sr = new SummarizeReplica();

            List<String> types = new ArrayList<>();
            List<Boolean> booleans = new ArrayList<>();
            booleans.add(true);
            booleans.add(false);
            types.add("tree");
            types.add("net");
            types.add("AZ");
            if(method.toLowerCase().equals("mcmc")) {
                sr.set_title_mcmc();
                for (String type : types) {
                    for (boolean bm : booleans) {
                        for (boolean bh : booleans) {
                            sr.summarize_mcmc(dir, type, bm, bh, burnin, sf);
                        }
                    }
                }
            }
            else if(method.toLowerCase().equals("ml")){
                sr.set_title_ML();
                List<String> gtt_list = new ArrayList<>();
                gtt_list.add("truegt");
                gtt_list.add("iqtree");

                for (String type : types) {
                    for (boolean bh : booleans) {
                        for (String gtt : gtt_list) {
                            sr.summarize_ML(dir, type, gtt, bh);
                        }
                    }
                }
            }

            writer.write(sr.sb.toString());

        } catch (FileNotFoundException e) {
            System.out.println(e.getMessage());
        }
    }





    public static void main(String[] args) {

        String output_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/replica2/mcmc_result.csv";
        String method = "MCMC";
//        String output_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/replica2/ML_result_dist.csv";
        int burnin = 20000000;
        int sf = 5000;
        checkTopologies(output_path, method, burnin, sf);

    }
}
