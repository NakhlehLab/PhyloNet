package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.test;
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
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SummarizeReplica {
//    Network truenet;
    StringBuilder sb = new StringBuilder();
    Map<String, Network> truenetmap = new HashMap<>();
    public static int chain_length = 20000000;
    public static int max_run_num = 3;

    /* Constructor */
    public SummarizeReplica(String method) {
//        truenetmap.put("tree", Networks.readNetwork("(G:5.0,(L:2.0,(Q:1.0,A:1.0):1.0):3.0);"));
//        truenetmap.put("net", Networks.readNetwork("(((Q:0.5)I8#H1:3.5::0.3,(L:2.0,(I8#H1:0.5::0.7,A:1.0)I4:1.0)I3:2.0)I2:1.0,G:5.0)I0;"));
//        truenetmap.put("AZ", Networks.readNetwork("((((A:0.1,B:0.1):0.2,C:0.3):0.1,D:0.4):10,E:10.4);"));
        if (method.toLowerCase().equals("ml")){
            truenetmap.put("net", Networks.readNetwork("(((((Q:0.75,#H1:0.25::0.3):0.75,R:1.5):2.5,(L:0.5)#H1:3.5::0.7):1,C:5),Z);"));
            truenetmap.put("tree", Networks.readNetwork("((((Q:1.5,R:1.5):2.5,L:4):1,C:5),Z);"));

        }
        else{
            truenetmap.put("net", Networks.readNetwork("((((Q:0.75,#H1:0.25::0.3):0.75,R:1.5):2.5,(L:0.5)#H1:3.5::0.7):1,C:5);"));
            truenetmap.put("tree", Networks.readNetwork("(((Q:1.5,R:1.5):2.5,L:4):1,C:5);"));
        }




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

    public Set<String> listFilesUsingJavaIO(String dir) {
        return Stream.of(new File(dir).listFiles())
                .filter(file -> !file.isDirectory())
                .map(File::getName)
                .collect(Collectors.toSet());
    }

    public String getLastLine(String path){
        String lastLine = "";
        try{

            BufferedReader br = new BufferedReader(new FileReader(path));
            String sCurrentLine = br.readLine();
            while (sCurrentLine != null)
            {
                lastLine = sCurrentLine;
                sCurrentLine = br.readLine();
            }
        }catch (Exception e){

        }
        return lastLine;
    }

    public void summarize_mcmc(String dir, String type, boolean murate, boolean homo, int locus_length, int burnin, int sf, String sum_type){
        Network true_net = truenetmap.get(type);
        int correct_cnt = 0;

        for (int i = 1; i <= 10; i++){
            String path = dir+"/"+i+"/";
            if (homo){
                path += "homo/";
            }
            else{
                path += "heter/";
            }
            path += locus_length +"/";
            if (murate){
                path += "murate/";
            }
            List<String> mcmc_path_list = new ArrayList<>();
//            String mcmc_path = "";

            for (int j = mcmc_path_list.size(); j < max_run_num; j++){
                for (String file : listFilesUsingJavaIO(path)){
                    if (file.startsWith("slurm-")){
                        String lastline = getLastLine(path+file);
                        if (lastline.contains(String.valueOf(chain_length)+"_"+String.valueOf(j))){
                            mcmc_path_list.add(path + file);
                        }
                    }
                }
            }

//            path += "mcmc.out";

            checkResultMCMC cr = new checkResultMCMC(burnin, sf);
            int cnt = 0;
            for(String mcmc_path : mcmc_path_list){
                cnt += cr.read_mcmc_out(mcmc_path);
            }

//            System.out.println(mcmc_path);
            Network net = null;
            Map<Network,Integer> net_list = cr.get_network_list();
            Map.Entry<Network, Integer> e = net_list.entrySet().iterator().next();
            if (sum_type.equals( "MAP")){
                net = cr.get_MAP_net();
            }
            else{

                net = e.getKey();
            }

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

    public void summarize_ML(String dir, String type, int locus_length, String gtt, boolean homo){
        if(!homo && gtt.equals("true")) {
            return;
        }
        Network true_net = truenetmap.get(type);

        for (int i = 1; i <= 10; i++){
            String path = dir+"/"+i+"/";
            if (homo){
                path += "homo/";
            }
            else{
                path += "heter/";
            }
            path += locus_length +"/ML_"+gtt+"_";



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
                sb.append(gtt.equals("true"));
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


    public static void  checkTopologies(String dir_path, String output_path, String method, int burnin, int sf, String net_type, String sum_type){
        try (PrintWriter writer = new PrintWriter(new File(output_path))) {
//            String dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/replica2/";
            SummarizeReplica sr = new SummarizeReplica(method);

            List<String> types = new ArrayList<>();
            List<Boolean> booleans = new ArrayList<>();
            booleans.add(true);
            booleans.add(false);
//            types.add("tree");
            types.add(net_type);
            int locus_length =2000;

//            types.add("AZ");
            if(method.toLowerCase().equals("mcmc")) {
                sr.set_title_mcmc();
                for (String type : types) {
                    for (boolean bm : booleans) {
                        for (boolean bh : booleans) {
                            sr.summarize_mcmc(dir_path, type, bm, bh, locus_length, burnin, sf, sum_type);
                        }
                    }
                }
            }
            else if(method.toLowerCase().equals("ml")){
                sr.set_title_ML();
                List<String> gtt_list = new ArrayList<>();
                gtt_list.add("true");
                gtt_list.add("iq");

                for (String type : types) {
                    for (boolean bh : booleans) {
                        for (String gtt : gtt_list) {
                            sr.summarize_ML(dir_path, type, locus_length, gtt, bh);
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
//        String dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1_2/medium/";
        String dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/net1_2/medium/";
        String output_path = dir_path+"/ML_result.csv";
        String sum_type = "MAP";
//        String output_path = dir_path+"/mcmc_result"+sum_type+".csv";
//        String method = "mcmc";
        String method = "ml";
//        String output_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/replica2/ML_result_dist.csv";
        int burnin = 2000000;
        int sf = 5000;

        checkTopologies(dir_path, output_path, method, burnin, sf, "net", sum_type);

    }
}
