package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by jz55 on 02/22/17.
 */
abstract class CommandBaseFileOutMatrix extends CommandBaseFileOut {

    protected Map<String, String> sourceIdentToMatrixData;

    public CommandBaseFileOutMatrix(SyntaxCommand motivatingCommand,
                                    ArrayList<Parameter> params,
                                    Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                    Proc3<String, Integer, Integer> errorDetected,
                                    RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    protected void parseMatrixData(String file) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            sourceIdentToMatrixData = new HashMap<>();
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

            while((s = br.readLine().trim()) != null) {
                if(s.toLowerCase().startsWith("format")) continue;
                if(s.toLowerCase().startsWith("matrix")) continue;
                if(s.toLowerCase().endsWith("end;")) {
                    break;
                }
                ss = s.trim().split("\\s+");
                if(ss.length == 2) {
                    sourceIdentToMatrixData.put(ss[0], ss[1]);
                }
            }
            if(sourceIdentToMatrixData.size() == 0) {
                sourceIdentToMatrixData = null;
            }
            br.close();

            if(sourceIdentToMatrixData.size() != nTaxon ) {
                throw new RuntimeException("\nWrong number of taxa\nIndicated: " + nTaxon + "\nRead: " + sourceIdentToMatrixData.size() + "\n");
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            sourceIdentToMatrixData = null;
        }
    }

    protected boolean assertDataExists(String taxaName, int line, int col) {
        if(!sourceIdentToMatrixData.containsKey(taxaName)) {
            errorDetected.execute(String.format("Unknown identifier '%s'.", taxaName), line, col);
            return false;
        }
        return true;
    }

}
