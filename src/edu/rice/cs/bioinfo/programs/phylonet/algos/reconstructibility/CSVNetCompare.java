package edu.rice.cs.bioinfo.programs.phylonet.algos.reconstructibility;
/*
 * @ClassName:   CSVNetCompare
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        10/6/23 3:20 PM
 */

import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.util.*;

public class CSVNetCompare {
    /* Constructor */
    private Map<String, Network> _netMap;

    /* Constructor */
    public CSVNetCompare(Map<String, Network> netMap) {
        _netMap = netMap;
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
            System.out.println(cur_net);
            System.out.println(value0_arr[1]);
            System.out.println(_netMap.get(value0_arr[1]));
            Network true_network = _netMap.get(value0_arr[1]);
            res.add(String.valueOf(Networks.computeDistanceBetweenTwoNetworks(cur_net, true_network)));
            Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(cur_net, true_network);
            res.add(String.valueOf(closest.Item3));
            res.add(String.valueOf(true_network.getReticulationCount()));
            res.add(String.valueOf(cur_net.getReticulationCount()));
            res.add(String.valueOf(closest.Item2.getReticulationCount()));
            res.add(String.valueOf(closest.Item1.getReticulationCount()));


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
            writer.write("difficulty,scale,num_gt,num_alleles,replica,max_reti_set,net,likelihood,distance,backbone_distance,true_net_reti,inferred_net_reti,true_backbone_reti,inferred_backbone_reti");
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

    public static void run(){
        String[] level_list = {"net_easy", "net_mid", "net_hard"};

//            for (String scale :scale_list){
//                for (int numgt: num_gts){
//                    for (int allele: num_alleles) {

        String dir_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/reconstructibility/data/";
        String method  = "MPL";
        String path = dir_path+ method+"_result.csv";

        String out_path = dir_path + method+"_dist.csv";
        String true_net = "";

        Map<String, Network> true_net_map = new HashMap<>();

        true_net = "(((A:2.5,((B:1,(C:0.5)#H1:0.5::0.7):1)#H2:0.5::0.6):1.5,(#H1:1::0.3,D:1.5):2.5):1,(#H2:1::0.4,E:3):2);";
        true_net_map.put("net_easy", Networks.readNetwork(true_net));

        true_net = "(((A:2.5,((B:1,(C:0.5)#H1:0.5::0.7):1)#H2:0.5::0.6):1.5,((#H1:1::0.3,D:1.5):1.5,#H2:1::0.4):1):1,E:5);";
        true_net_map.put("net_mid", Networks.readNetwork(true_net));

        true_net = "(((A:2.5,(B:1,((C:0.5)#H1:0.25::0.7)#H2:0.25::0.6):1.5):1.5,((#H1:1::0.3,D:1.5):1.5,#H2:2.25::0.4):1):1,E:5);";
        true_net_map.put("net_hard", Networks.readNetwork(true_net));

        CSVNetCompare cn = new CSVNetCompare(true_net_map);
        List<List<String>> records = cn.read_csv(path);
        cn.write_csv(records, out_path);


    }

    public static void main(String[] args) {
        run();
    }

}
