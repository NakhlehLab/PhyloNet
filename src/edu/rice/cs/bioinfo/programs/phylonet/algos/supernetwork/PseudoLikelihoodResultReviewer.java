package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import com.google.gson.Gson;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 1/22/19
 * Time: 2:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class PseudoLikelihoodResultReviewer {
    public static void main(String args[]) {
        String networkfilename = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/Networks/networks.json";
        SimTest.NetworkListJson networkListjson = null;
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(networkfilename));

            Gson gson = new Gson();
            networkListjson = gson.fromJson(bufferedReader, SimTest.NetworkListJson.class);

        } catch(Exception e) {
            e.printStackTrace();
        }

        List<String> netstrings = new ArrayList<>();
        for(SimTest.NetworkListJson.NetworkJson networkJson : networkListjson.networks) {
            System.out.println(networkJson.tag);
            netstrings.add(networkJson.netstring);
        }


        //String resultPrefix = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run37/slurm-6159678_";
        String resultPrefix = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run68/slurm-7114606_";

        for(int i = 0; i < 24 ; i++) {
            System.out.println(i);

            String filename = resultPrefix + i + ".out";
            Network inferredNetwork = null;
            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                while((s = in.readLine()) != null) {
                    if(s.equals("Inferred Network #1:")) {
                        s = in.readLine();
                        inferredNetwork = Networks.readNetwork(s);
                        break;
                    }

//                    if(s.equals("Inferred Network #1:")) {
//                        s = in.readLine();
//                        inferredNetwork = Networks.readNetwork(s);
//                        break;
//                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

            Network trueNetwork = Networks.readNetwork(netstrings.get(i));
            System.out.println("True: " + Networks.getDendroscopeCompatibleString(trueNetwork));
            System.out.println("Inferred: " + inferredNetwork);
            System.out.println("Inferred # Reti: " + inferredNetwork.getReticulationCount());

            Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(inferredNetwork, trueNetwork);
            System.out.println("Closest true # Reti: " + closest.Item2.getReticulationCount());
            System.out.println(closest.Item2);
            System.out.println("Closest inferred # Reti: " + closest.Item1.getReticulationCount());
            System.out.println(closest.Item1);
            System.out.println("Distance: " + closest.Item3);

            System.out.println();
        }
    }
}
