package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 3/23/18
 * Time: 4:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimSeqInGTBySeqGen {
    static private String getFreqString(double freq[]) {
        return " -f" + freq[0] + "," + freq[1] + "," + freq[2] + "," + freq[3] + " ";
    }

    static private String getRateString(double rates[]) {
        return " -r" + rates[0] + "," + rates[1] + "," + rates[2] + "," + rates[3] + "," + rates[4] + "," + rates[5] + " ";
    }

    static private String getLenString(int len) {
        return " -l" + len + " ";
    }

    static private String getSeedString(int seed) {
        return " -z" + seed + " ";
    }

    static public Map<String, String> execute(Tree gt, double theta, double freq[], double rates[], int len, Integer seed, String SeqGenPath) {
        StringBuilder sb = new StringBuilder();
        sb.append(SeqGenPath);

        sb.append(" -mGTR ");
        sb.append("-s");
        sb.append(theta / 2.0);
        sb.append(getFreqString(freq));
        sb.append(getRateString(rates));
        sb.append(getLenString(len));
        if(seed != null) {
            sb.append(getSeedString(seed));
        }

        try {
            PrintWriter out = new PrintWriter("seqgen.tree");
            out.print(gt.toNewick());
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        sb.append(" seqgen.tree ");

        String command = sb.toString();
        Map<String, String> seq = new HashMap<>();

        try{
            Process proc = Runtime.getRuntime().exec(command,null,null);

            BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            String line;
            int retVal = proc.waitFor();

            boolean skip = true;
            while((line=stdout.readLine())!=null){
                line = line.trim();
                if(line.length() == 0) {
                    continue;
                }

                if(skip) {
                    skip = false;
                    continue;
                }

                String parts[] = line.split(" ");
                seq.put(parts[0], parts[parts.length - 1]);
            }
            stdout.close();
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }

        return seq;

    }
}
