package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.test;
/*
 * @ClassName:   checkResultML
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        10/13/21 7:26 PM
 */

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

public class checkResultML {
    private Network topology = null;
    private double likelihood = 0.0;

    /* Constructor */
    public checkResultML(String path) {
        if (!(new File(path)).exists()){
            return;
        }
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String s;

            boolean begin = false;

            while((s = br.readLine()) != null) {
                s = s.trim();
                if(s.contains("Inferred Network #1:")) {
                    begin = true;
                    continue;
                }

                if(begin) {
                    if (s.startsWith("(")){
                        topology = Networks.readNetwork(s);
                    }
                    else if(s.startsWith("Total log probability:")){
                        likelihood = Double.parseDouble(s.split(": ")[1].trim());
                        break;
                    }

                }

            }

            br.close();

        } catch(Exception e) {
            e.printStackTrace();
        }

    }

    public Network getTopology(){
        return topology;

    }

    public double getLikelihood(){
        return likelihood;
    }





    public static void main(String[] args) {
        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/tree_sim/experiment/butterfly/mcmc/original/varyps/mcmc.out";
//        checkResultMCMC cr = new checkResultMCMC(path);
    }
}
