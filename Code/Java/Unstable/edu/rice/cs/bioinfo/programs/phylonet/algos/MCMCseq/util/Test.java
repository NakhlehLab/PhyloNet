package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 5/5/18
 * Time: 8:57 AM
 * To change this template use File | Settings | File Templates.
 */
public class Test {
    static void convertMSTrees(String filename, double scale) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            PrintWriter out = new PrintWriter(filename + ".cvt");
            String line;
            int currentNumber = 1;
            String numbers = "";
            while ((line = br.readLine()) != null) {
                Network net = Networks.readNetwork(line);
                for(Object nodeObj : Networks.postTraversal(net)) {
                    NetNode node = (NetNode) nodeObj;
                    for(Object parentObj : node.getParents()) {
                        NetNode parent = (NetNode) parentObj;
                        node.setParentDistance(parent, node.getParentDistance(parent) * scale);
                    }
                }
                String treeName = "gt" + currentNumber;
                if(currentNumber == 1) {
                    numbers = treeName;
                } else {
                    numbers = numbers + "," + treeName;
                }
                out.println("Tree " + treeName + " = " + net);
                currentNumber++;
            }
            out.close();

            System.out.println("(" + numbers + ")");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        convertMSTrees("/Users/zhujiafan/Documents/BioinfoData/EZONChapter/genetrees_0.txt", 0.01);
    }
}
