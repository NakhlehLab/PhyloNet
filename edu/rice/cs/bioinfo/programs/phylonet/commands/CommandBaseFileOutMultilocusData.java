package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

/**
 * Created by dw20 on 10/27/16.
 */
abstract class CommandBaseFileOutMultilocusData extends CommandBaseFileOut {

    protected Map<String, Map<String, String>> sourceIdentToMultilocusData;

    public CommandBaseFileOutMultilocusData(SyntaxCommand motivatingCommand,
                                    ArrayList<Parameter> params,
                                    Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                    Proc3<String, Integer, Integer> errorDetected,
                                    RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    protected void parseMultiLociData(String file) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            sourceIdentToMultilocusData = new HashMap<>();
            String s;
            String[] ss;
            int nTaxon = 0;
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
                    nTaxon = Integer.parseInt(s.substring(start, end));
                }
                if(begin && nTaxon > 0) break;
            }

            s = br.readLine().trim();
            while(s != null) {
                if(s.toLowerCase().endsWith("end;")) {
                    break;
                }
                if(s.startsWith("[") && s.endsWith("]")) {
                    ss = s.substring(1, s.length() - 1).replaceAll("\\s+", "").split(",");
                    if(ss.length < 2) {
                        br.close();
                        throw new RuntimeException(s);
                    }
                    String name = ss[0];
                    int length = Integer.parseInt(ss[1]);
                    Map<String, String> taxonSeqMap = new HashMap<>();
                    s = br.readLine().trim();
                    while(!(s.startsWith("[") && s.endsWith("]")) && (ss = s.split("\\s+")).length == 2) {
                        if(ss[1].length() != length) {
                            throw new RuntimeException("wrong sequence length " + ss[1].length() + "\n" + s);
                        }
                        taxonSeqMap.put(ss[0], ss[1]);
                        s = br.readLine().trim();
                    }
                    sourceIdentToMultilocusData.put(name, taxonSeqMap);
                } else {
                    s = br.readLine().trim();
                }
            }
            if(sourceIdentToMultilocusData.size() == 0) {
                sourceIdentToMultilocusData = null;
            }
            br.close();
        } catch (Exception ex) {
            sourceIdentToMultilocusData = null;
        }
    }

    protected boolean assertDataExists(String lociName, int line, int col) {
        if(!sourceIdentToMultilocusData.containsKey(lociName)) {
            errorDetected.execute(String.format("Unknown identifier '%s'.", lociName), line, col);
            return false;
        }
        return true;
    }

}
