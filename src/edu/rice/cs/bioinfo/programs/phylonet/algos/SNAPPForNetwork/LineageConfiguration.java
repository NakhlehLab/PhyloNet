package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 3/7/19
 * Time: 2:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class LineageConfiguration {
    private Map<String, Map<String, R>> top;
    private Map<String, Map<String, R>> bottom;

    public LineageConfiguration() {
        top = new HashMap<>();
        bottom = new HashMap<>();
    }

    public void setTop(NetNode child, NetNode parent, R r) {
        if(!top.containsKey(child.getName())) top.put(child.getName(), new HashMap<>());
        top.get(child.getName()).put(parent.getName(), r);
    }

    public void setBottom(NetNode child, NetNode parent, R r) {
        if(!bottom.containsKey(child.getName())) bottom.put(child.getName(), new HashMap<>());
        bottom.get(child.getName()).put(parent.getName(), r);
    }

    public R getTop(NetNode child, NetNode parent) {
        return top.get(child.getName()).get(parent.getName());
    }

    public R getBottom(NetNode child, NetNode parent) {
        return bottom.get(child.getName()).get(parent.getName());
    }
}
