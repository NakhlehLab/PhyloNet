package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by zhiyan on Oct, 2020
 */
public class Tools {
    public static List<Alignment> readSeq(Map<String, Map<String, String>> multiLociData) {
        ArrayList<Alignment> alns = new ArrayList<>();
        if(Utils._TAXON_MAP != null){
            Set<String> alleleSet = new HashSet<>();
            for(String species :Utils._TAXON_MAP.keySet() ) {
                alleleSet.addAll(Utils._TAXON_MAP.get(species));
            }
            for (String locusName: multiLociData.keySet()){
                Map<String, String> locus = multiLociData.get(locusName);
                Map<String, String> newLocus = new HashMap<>();
                for(String alleleName : locus.keySet()) {
                    if(alleleSet.contains(alleleName)) {
                        newLocus.put(alleleName, locus.get(alleleName));
                    }
                }
                alns.add(new Alignment(newLocus, locusName));
            }
        }
        else{
            for (String locusName : multiLociData.keySet()) {
                alns.add(new Alignment(multiLociData.get(locusName), locusName));
            }
        }
        return alns;
    }

    public static Map<String, Map<String, String>> parseNexusFile(String file) {
        Map<String, Map<String, String>> multiLociData = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String s;
            String[] ss;
            int taxaNum = 0;
            boolean begin = false;
            while((s = br.readLine().trim()) != null ) {
                s = s.replaceAll("\\s+", "").toLowerCase();
                if (s.startsWith("begindata;")) {
                    begin = true;
                }
                if (s.contains("ntax=")) {
                    int start = s.indexOf("ntax=") + 5;
                    int end = start + 1;
                    while(end < s.length() && Character.isDigit(s.charAt(end))) {
                        end++;
                    }
                    taxaNum = Integer.parseInt(s.substring(start, end));
                }
                if(begin && taxaNum > 0) break;
            }

            s = br.readLine().trim();

            /*   read alignment   */
            while(s != null) {
                if (s.toLowerCase().endsWith("end;")) break;
                if (s.startsWith("[") && s.endsWith("]")) {
                    ss = s.substring(1, s.length() - 1).replaceAll("\\s+", "").split(",");
                    if(ss.length < 2) {
                        br.close();
                        throw new RuntimeException(s);
                    }
                    String locusName = ss[0];
                    int seqLength = Integer.parseInt(ss[1]);
                    Map<String, String> taxonSeqMap = new HashMap<>();
                    s = br.readLine().trim();
                    while(!(s.startsWith("[") && s.endsWith("]")) && (ss = s.split("\\s+")).length == 2) {
                        if(ss[1].length() != seqLength) {
                            throw new RuntimeException("wrong sequence length " + ss[1].length() + "\n" + s);
                        }
                        taxonSeqMap.put(ss[0], ss[1]);
                        s = br.readLine().trim();
                    }
                    multiLociData.put(locusName, taxonSeqMap);
                } else {
                    s = br.readLine().trim();
                }
            }
            if (multiLociData.size() == 0) { multiLociData = null; }


            /*   read command */

            s = br.readLine();
            while(s != null){
                s = s.trim();
//                if (s.toLowerCase().startsWith("begin")) continue;
                if (s.toLowerCase().startsWith("mcmc")) {
                    /* taxon map */
                    if(s.contains(" -tm ")) {
                        String sb = s.substring(s.indexOf('<')+1, s.indexOf('>'));
                        if (Utils._TAXON_MAP == null || Utils._TAXON_MAP.isEmpty()){
                            Utils._TAXON_MAP = new HashMap<>();
                        }
                        ss = sb.split(";");
                        if (sb.equals("") || sb == null) break;
                        for (String allelemapping :ss){
                            String[] arr = allelemapping.trim().split(":");
                            if (arr.length != 2) throw new IOException("wrong input format of taxonmap");

                            Utils._TAXON_MAP.put(arr[0], new ArrayList<>());
                            String[] alleles = arr[1].split(",");
                            for (int i = 0; i < alleles.length; i++){
                                Utils._TAXON_MAP.get(arr[0]).add(alleles[i]);
                            }
                        }
                    }
                    /* MC3 setting */

                    if(s.contains(" -mc3 ")){
                        String mc3para = s.substring(s.indexOf("mc3"));
                        String sb = mc3para.substring(mc3para.indexOf('(')+1, mc3para.indexOf(')'));
                        String [] arr = sb.split(",");
                        for (String temperature: arr){
                            Utils._MC3_CHAINS.add(Double.parseDouble(temperature.trim()));
                        }
                    }

                    if (s.contains("-varyps")){
                        Utils._ESTIMATE_POP_SIZE = true;
                        Utils._CONST_POP_SIZE = false;
                    }
                    else if (s.contains("-fixps")) {
                        Utils._ESTIMATE_POP_SIZE = false;
                        Utils._CONST_POP_SIZE = true;
                        String sb = s.substring(s.indexOf("-fixps")+6).trim();
                        sb = sb.split(" ")[0].trim();
                        if (sb.endsWith(";")){
                            sb = sb.substring(0, sb.length()-1);
                        }
                        System.out.println(sb);
                        Utils._POP_SIZE_MEAN = Double.parseDouble(sb.split(" ")[0].trim());
                    }

                    if (s.contains(" -pl ")){
                        System.out.println(s.indexOf(" -pl "));
                        String sb = s.substring(s.indexOf("-pl ")+3).trim();
                        sb = sb.split(" ")[0].trim();
                        Utils._NUM_THREADS = Integer.parseInt(sb);
                    }


                }
                s = br.readLine();
            }

            br.close();

        } catch (Exception e) {
            multiLociData = null;
        }
        return multiLociData;
    }

    public static long getSeed(int seedIndex) {
        long[] seedList = {
                1662750824,
                889504735,
                1469796979,
                570425800,
                2110580495,
                165088863,
                5546266,
                748834402,
                495282039,
                2050429893,
                1716904794,
                533376665,
                676979682,
                429867957,
                1150117749,
                1410492399,
                1767955838,
                570829003,
                1064459265,
                1184309624,
                1848926225,
                1197043800,
                1524544731,
                983013182,
                1387597549,
                1893716468,
                1214831113,
                1258906609,
                719783379,
                1242857737,
                397642604,
                1021577796,
                1124883814,
                1040570479,
                1419532573,
                2115293527,
                1742809030,
                1269281837,
                1016045416,
                1513715595,
                1804408343,
                492979769,
                1408038837,
                493339225,
                1573050344,
                936211922,
                2127368446,
                1381344219,
                1413941022,
                1698673633,
                343016179,
                1786043843,
                524271889,
                1164018220,
                337503320,
                574387154,
                1322290499,
                1127005390,
                1301322448,
                1117501455,
                139344831,
                1455293947,
                1731258013,
                1855578366,
                1689723492,
                1531382810,
                522955034,
                335131325,
                1715362749,
                1708931807,
                835282185,
                644232995,
                1723025642,
                14826937,
                121229473,
                1497059636,
                1400292868,
                1407865436,
                902037916,
                803351961,
                1291422308,
                2028656849,
                823281131,
                502106042,
                588062233,
                2043799255,
                1617895567,
                1925681424,
                870100593,
                966999316,
                348106262,
                903218539,
                939955831,
                647010041,
                575416952,
                1517578247,
                932200576,
                1831137971,
                1074748807,
                1474811250};
        return seedList[seedIndex];
    }

}
