package edu.rice.cs.bioinfo.programs.phylonet.algos.chi2;
/*
 * @ClassName:   Check_net
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        8/3/22 2:29 PM
 */
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.nio.file.Files;
import java.util.*;


public class Check_net {
    private Network _net;

    /* Constructor */
    public Check_net(String true_net) {
        _net = Networks.readNetwork(true_net);
    }

    private List<String> getRecordFromLine(String line) {
        List<String> values = new ArrayList<String>();
        List<String> res = new ArrayList<String>();

        try (Scanner rowScanner = new Scanner(line)) {
            rowScanner.useDelimiter("\"");
            while (rowScanner.hasNext()) {
                String x = rowScanner.next();
                values.add(x);
            }
        }
//        System.out.println(values.size());
        if (values.size() == 3){
            String[] value0_arr = values.get(0).split(",");
            res.add(value0_arr[1]);
            res.add(value0_arr[2]);
            res.add(value0_arr[3]);
            res.add(value0_arr[4]);
            res.add(value0_arr[5]);
            res.add(value0_arr[6]);
//            System.out.println(values.get(1));
            String net_str = values.get(1);
            net_str = net_str.replaceAll("\\[.*?\\]","");
//            System.out.println(net_str);
            Network cur_net = Networks.readNetwork(net_str);
            res.add("\""+net_str+"\"");

            String[] arr = values.get(2).split(",");
            res.add(arr[1]);
//            res.add(arr[2]);
            res.add(String.valueOf(Networks.computeDistanceBetweenTwoNetworks(cur_net, _net)));

        }
//        System.out.println(res.size());
        return res;
    }

    public List<List<String>> read_csv(String csv_path){
        List<List<String>> records = new ArrayList<>();
        try (Scanner scanner = new Scanner(new File(csv_path));) {
            while (scanner.hasNextLine()) {
                records.add(getRecordFromLine(scanner.nextLine()));
            }
            System.out.println(records);


        }catch (FileNotFoundException e) {
            System.out.println(e.getMessage());
        }
        return records;
    }

    public void write_csv(List<List<String>> records, String out_path){
        try{
            BufferedWriter writer = new BufferedWriter(new FileWriter(out_path));
            writer.write("difficulty,scale,num_gt,num_alleles,replica,max_reti_set,net,likelihood,distance");
            for (List<String> line : records){
                writer.write(String.join(",", line));
                writer.write("\n");
            }

            writer.flush();
            writer.close();
        }catch (IOException e){
            System.out.println(e.getMessage());
        }
    }

    public static void main(String[] args) {
//        String[] level_list = {"net0", "net1"};
        String[] level_list = {"net1_2"};
        String[] scale_list = {"long"};
        String[] heter_list = {"homo", "heter"};
        String[] murate_list = {"_murate", ""};
        for(String level : level_list){
            for (String scale :scale_list){
                for (String heter: heter_list){
                    for (String murate: murate_list){
                        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/"+level+"/"+scale+"/posterior_prob_trace_"+heter+murate+".csv";
                        String out_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mcmc/"+level+"/"+scale+"/posterior_prob_trace_"+heter+murate+"_out.csv";
                        String true_net = "";
                        if (level.equals("net0")){
                            true_net = "(((Q:0.6,R:0.6)I3:1.0,L:1.6)I1:0.4,C:2);";
                        }
                        else if (level.equals("net1")){
                            true_net = "((((Q:0.4)#H1:0.2::0.7,R:0.6)I3:1.0,(L:0.8,#H1:0.4::0.3):0.8)I1:0.4,C:2);";
                        }
                        else if (level.equals("net1_2")){
                            true_net = "((((Q:0.3,#H1:0.1::0.3):0.3,R:0.6):1,(L:0.2)#H1:1.4::0.7):0.4,C:2);";
                        }
                        else if (level.equals("net1_3")){
                            true_net = "((((Q:0.3,#H1:0.1::0.7):0.3,R:0.6):1,(L:0.2)#H1:1.4::0.3):0.4,C:2);";
                        }
                        else if (level.equals("net1_4")){
                            true_net = "((((Q:0.3,#H1:0.1::0.7):0.3,R:0.6):1,(L:0.2)#H1:1.4::0.3):0.4,C:2);";
                        }

                        Check_net cn = new Check_net(true_net);
                        List<List<String>> records = cn.read_csv(path);
                        cn.write_csv(records, out_path);
                    }

                }
            }
        }


    }
}
