package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import com.google.gson.Gson;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 1/22/19
 * Time: 10:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class RunningTimeReviewer {
    static Map<String, Tuple<List<String>, Double>> visitAllRunningTime(String folder) {
        File path = new File(folder);

        if(!path.exists()) return null;

        List<String> filenames = new ArrayList<>();
        Map<String, Tuple<List<String>, Double>> result = new HashMap<>();

        File [] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].toString().endsWith(".out")) {
                    //System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);

        for(String filename : filenames) {
            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                double time = 3600 * 4;
                int total = 0;
                int end = 400;
                boolean begin = false;
                Network curSample = null;

                while((s = in.readLine()) != null) {
                    if (begin) {
                        if (s.startsWith("[")) {
                            curSample = Networks.readNetworkWithRootPop(s);
                            total++;

                        } else if(s.contains("Summarization")) {
                            begin = false;
                        }
                    } else {
                        if (s.contains("Logger")) {
                            begin = true;
                        } else if (s.startsWith("Total elapsed time : ")) {
                            time = Double.parseDouble(s.substring(s.indexOf(":") + 2, s.length() - 2));
                            break;
                        }
                    }
                }

                if(curSample == null) continue;

                List<String> leafname = new ArrayList<>();
                for(Object leafObj : curSample.getLeaves()) {
                    NetNode leaf = (NetNode) leafObj;
                    leafname.add(leaf.getName());
                }
                Collections.sort(leafname);

                if(total != end) {
                    time = time / total * end;
                    //System.out.println(filename + " " + total);
                }

                result.put(filename, new Tuple<>(leafname, time));


            } catch (IOException e) {
                e.printStackTrace();
                return null;
            }
        }
        return result;
    }

    static void checkLizardTime() {
        String folder = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run72";
        File path = new File(folder);

        List<String> filenames = new ArrayList<>();
        Map<String, Tuple<List<String>, Double>> result = new HashMap<>();

        File [] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].toString().endsWith(".out")) {
                    //System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);

        for(String filename : filenames) {
            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                double time = 3600 * 24;
                int total = 0;
                int end = 1200;
                boolean begin = false;
                Network curSample = null;

                while((s = in.readLine()) != null) {
                    if (begin) {
                        if (s.startsWith("[")) {
                            curSample = Networks.readNetworkWithRootPop(s);
                            total++;

                        } else if(s.contains("Summarization")) {
                            begin = false;
                        }
                    } else {
                        if (s.contains("Logger")) {
                            begin = true;
                        } else if (s.startsWith("Total elapsed time : ")) {
                            time = Double.parseDouble(s.substring(s.indexOf(":") + 2, s.length() - 2));
                            break;
                        }
                    }
                }

                if(curSample == null) continue;

                List<String> leafname = new ArrayList<>();
                for(Object leafObj : curSample.getLeaves()) {
                    NetNode leaf = (NetNode) leafObj;
                    leafname.add(leaf.getName());
                }
                Collections.sort(leafname);

                if(total != end) {
                    time = time / total * end;
                    //System.out.println(filename + " " + total);
                }

                result.put(filename, new Tuple<>(leafname, time));


            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        double total = 0;
        for(String filename : result.keySet()) {
            total += result.get(filename).Item2;
        }
        System.out.println("Total time (s): " + total);
        System.out.println("Total CPU-Hours: " + total / 3600 * 2); // 2 cores
    }

    public void ReduceTrinets(Map<String, String> allele2species, String list_path) {
        Set<Set<String>> triset = new HashSet<>();
        String outgroup = "Z";
        try {
            BufferedReader in = new BufferedReader(new FileReader(list_path));
            String s;
            while((s = in.readLine()) != null) {
                List<String> cur_alleles = new ArrayList<>();
                String cur_allele = "";
                Set<String> cur_species = new HashSet<>();
                for(int i = 0 ; i < s.length() ; i++) {
                    if(s.charAt(i) == '[' || s.charAt(i) == ' ') {
                        continue;
                    } else if(s.charAt(i) == ',' || s.charAt(i) == ']') {
                        cur_alleles.add(cur_allele);
                        cur_allele = "";
                    } else {
                        cur_allele += s.charAt(i);
                    }
                }

                for(String allele : cur_alleles) {
                    cur_species.add(allele2species.get(allele));
                }

                if(cur_species.size() == 1) {
                    continue;
                } else if(cur_species.size() == 2) {
                    cur_species.add(outgroup);
                }

                triset.add(cur_species);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }


//        subnetworks_reduceed_ = new ArrayList<>();
//        Iterator<NetworkWithInfo> it = subnetworks_.iterator();
//        while (it.hasNext()) {
//
//            NetworkWithInfo netinfo = it.next();
//
//            Set<String> taxa = new HashSet<>(netinfo.taxa);
//            if(!triset.contains(taxa) /*&& !taxa.contains(outgroup)*/) {
//                it.remove();
//                subnetworks_reduceed_.add(netinfo);
//            }
//        }
//
//        if(printDetails_) {
//            System.out.println("Reduce number of trinets to " + subnetworks_.size());
//        }



    }

    static Tuple<List<List<String>>, Double> visitReducedTrinets(String stage2file) {
        double time = 0;
        List<List<String>> trinets = new ArrayList<>();

        try {
            BufferedReader in = new BufferedReader(new FileReader(stage2file));
            String s;

            while((s = in.readLine()) != null) {
                if(s.startsWith("Time (s): ")) {
                    time = Double.parseDouble(s.substring(s.indexOf(": ") + 2));
                } else if(s.startsWith("After check, number of trinets is")) {

                }
            }

        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        return  null;
    }

    static Integer visitReducedSize(String stage2file) {
        double time = 0;
        List<List<String>> trinets = new ArrayList<>();
        int num = 0;
        int count = 0;

        try {
            BufferedReader in = new BufferedReader(new FileReader(stage2file));
            String s;

            while((s = in.readLine()) != null) {
                if(s.startsWith("Reduce number of trinets to")) {
                    count = 1;
                    num = Integer.parseInt(s.substring(s.lastIndexOf(" ") + 1));
                } else if(s.startsWith("Need more trinets: ")) {
                    count++;
                    num += Integer.parseInt(s.substring(s.lastIndexOf(" ") + 1));
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        System.out.println("\t" + num + "\t" + count);
        return  null;
    }

    static Double visitStage2Time(String stage2file) {
        double time = 0;
        List<List<String>> trinets = new ArrayList<>();

        try {
            BufferedReader in = new BufferedReader(new FileReader(stage2file));
            String s;

            while((s = in.readLine()) != null) {
                if(s.startsWith("Time (s): ")) {
                    return Double.parseDouble(s.substring(s.lastIndexOf(" ") + 1));
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        return  null;
    }

    public static void checkAllTime(String args[]) {
        String networkfilename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/networks.json";
        String resultRoot = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run61/"; // run61 run47

        SimTest.NetworkListJson networkListjson = null;
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(networkfilename));

            Gson gson = new Gson();
            networkListjson = gson.fromJson(bufferedReader, SimTest.NetworkListJson.class);

        } catch(Exception e) {
            e.printStackTrace();
        }

        int index = 0;
        int requiredIndex = -1;
        int numCores = 2;
        for(SimTest.NetworkListJson.NetworkJson networkJson : networkListjson.networks) {
            if(index++ < requiredIndex) {continue;}
            String resultPath = resultRoot + "/" + networkJson.tag;
            Map<String, Tuple<List<String>, Double>> info = visitAllRunningTime(resultPath);
            double allTime = 0.0;
            for(String filename : info.keySet()) {
                allTime += info.get(filename).Item2;
            }
            System.out.println(networkJson.tag + "\t" + allTime / 3600.0 * numCores); //cpu-hours

            //break;
        }
    }

    public static void checkAllReduced(String args[]) {
        String networkfilename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/networks.json";
        String resultRoot = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run61/"; // run61 run47
        String stage2root = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/simtest_revison3_reduce/";

        SimTest.NetworkListJson networkListjson = null;
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(networkfilename));

            Gson gson = new Gson();
            networkListjson = gson.fromJson(bufferedReader, SimTest.NetworkListJson.class);

        } catch(Exception e) {
            e.printStackTrace();
        }

        int index = -1;
        int requiredIndex = -1;
        for(SimTest.NetworkListJson.NetworkJson networkJson : networkListjson.networks) {
            if(index++ < requiredIndex) {continue;}
            String resultFile = stage2root + "/simtest_reduce_run61_" + index + ".txt";
            visitReducedSize(resultFile);

            //System.out.println(networkJson.tag + "\t" + visitStage2Time(resultFile));

            //break;
        }
    }

    public static void checkAllFull(String args[]) {
        String networkfilename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/networks.json";
        String resultRoot = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run61/";
        String stage2root = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/simtest_revison3/";

        SimTest.NetworkListJson networkListjson = null;
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(networkfilename));

            Gson gson = new Gson();
            networkListjson = gson.fromJson(bufferedReader, SimTest.NetworkListJson.class);

        } catch(Exception e) {
            e.printStackTrace();
        }

        int index = -1;
        int requiredIndex = -1;
        for(SimTest.NetworkListJson.NetworkJson networkJson : networkListjson.networks) {
            if(index++ < requiredIndex) {continue;}
            String resultFile = stage2root + "/simtest_run61_" + index + ".txt";

            System.out.println(networkJson.tag + "\t" + visitStage2Time(resultFile));

            //break;
        }
    }

    public static void main(String args[]) {
        checkAllReduced(args);
        //checkAllFull(args);
        //checkLizardTime();
    }
}
