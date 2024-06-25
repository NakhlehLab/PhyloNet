package edu.rice.cs.bioinfo.programs.phylonet.algos.chi2;
/*
 * @ClassName:   UltraNetChecker
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        3/22/22 10:57 PM
 */

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;

public class UltraNetChecker {
    /* Constructor */
    public UltraNetChecker() {

    }

    public static void main(String[] args) {
        String netstr = "(((((Q:0.2)#H1:0.1::0.7,R:0.3)I3:0.5,(L:0.4,#H1:0.3::0.3):0.4)I1:0.2,(G:0.4,C:0.4)I2:0.6)I0:1,Z:2);";
        UltrametricNetwork network = new UltrametricNetwork(netstr, new ArrayList<>());
    }
}
