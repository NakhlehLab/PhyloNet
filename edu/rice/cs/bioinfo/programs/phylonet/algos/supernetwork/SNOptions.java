package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/21/18
 * Time: 2:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class SNOptions {
    public double eps = 0.01;
    public boolean trustReticulationTime = true;
    public boolean reconcileHeights = false;
    public String outgroup = "Z";
    public Network trueNetwork = null;
    public String tripletFilename = null;
    public Map<String, String> allele2species = null;
    public List<String> buildOrder = null;
    public Set<String> backboneLeaves = null;
}
