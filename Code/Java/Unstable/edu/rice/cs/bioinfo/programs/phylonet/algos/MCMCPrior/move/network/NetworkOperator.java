package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * Created by wendingqiao on 2/15/16.
 * Abstract operator that changes network
 */
public abstract class NetworkOperator extends Operator {


    protected UltrametricNetwork _network;
    protected boolean _violate; // if the operation would violate the temporal constraint

    public NetworkOperator(UltrametricNetwork net)  {
        this._network = net;
        this._violate = true; // all moves will violate constraints except node height & pop size
    }

    public Utils.MOVE_TYPE getCategory() {
        return Utils.MOVE_TYPE.NETWORK;
    }

    public boolean mayViolate() {
        return  _violate;
    }

    public void optimize (double logAlpha) {};

    public void setNodeHeight(NetNode<NetNodeInfo> node, double height) {
        node.getData().setHeight(height);
        for(NetNode<NetNodeInfo> child : node.getChildren()) {
            child.setParentDistance(node, height - child.getData().getHeight());
        }
        for(NetNode<NetNodeInfo> par : node.getParents()) {
            node.setParentDistance(par, par.getData().getHeight() - height);
        }
    }

    protected void adopt(NetNode<NetNodeInfo> par, NetNode<NetNodeInfo> child, double[] params) {
        par.adoptChild(child, par.getData().getHeight() - child.getData().getHeight());
        child.setParentProbability(par, params[0]);
        child.setParentSupport(par, params[1]);
    }

    protected double[] getParameters(NetNode<NetNodeInfo> par, NetNode<NetNodeInfo> child) {
        return new double[] {child.getParentProbability(par), child.getParentSupport(par)};
    }

}
