package gridSearch;

import phylogeny.Node;

public class EvoTreeNob extends Nob {

    private Node childNode;

    public EvoTreeNob(int gIn, double minIn, double maxIn, Node childNodeIn)  {
        super (gIn, minIn, maxIn);
        this.childNode = childNodeIn;
    }

    @Override
    public void set_param(double value) {
        childNode.setTbranch(value);
    }

    @Override
    public double get_param() {
        return childNode.getTbranch();
    }

}
