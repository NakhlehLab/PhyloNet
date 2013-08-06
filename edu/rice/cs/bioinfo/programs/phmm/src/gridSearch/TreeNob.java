package gridSearch;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

public class TreeNob extends Nob {
	
	private TNode childNode;
	
	public TreeNob(int gIn, double minIn, double maxIn, TNode childNodeIn)  {
		super (gIn, minIn, maxIn);
		this.childNode = childNodeIn;
	}

	@Override
	public void set_param(double value) {
		childNode.setParentDistance(value);
	}

	@Override
	public double get_param() {
		return childNode.getParentDistance();
	}

}