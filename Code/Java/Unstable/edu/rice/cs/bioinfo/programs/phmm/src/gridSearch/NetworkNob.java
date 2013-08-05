package gridSearch;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

/**
 * This is a the knob used for the tree branch parameters from 
 * Network trees.
 * 
 * @author k3kathy
 *
 */
public class NetworkNob extends Nob {
	
	private NetNode<Double> netNode;
	private NetNode<Double> parentNode;
	
	public NetworkNob(int gIn, double minIn, double maxIn, NetNode<Double> netNodeIn, NetNode<Double> parentNodeIn) {
		super(gIn, minIn, maxIn);
		
		this.netNode = netNodeIn;
		this.parentNode = parentNodeIn;
	}
	
	public void set_param(double value) {
		netNode.setParentDistance(parentNode, value);
	}
	
	public double get_param() {
		return netNode.getParentDistance(parentNode);
	}
	
	
}
