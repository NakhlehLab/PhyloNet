package gridSearch;

//import phylogeny.Node;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

public class GeneGenealogyNob extends Nob {

    private TNode childNode;

    public GeneGenealogyNob (int gIn, double minIn, double maxIn, TNode childNodeIn)  {
        super (gIn, minIn, maxIn);
        this.childNode = childNodeIn;
    }

    @Override
    public void set_param(double value) {
        backupParam = childNode.getParentDistance();
        childNode.setParentDistance(value);
    }

    @Override
    public double get_param() {
        return childNode.getParentDistance();
    }

    public void restoreParameterValue() {
        childNode.setParentDistance(backupParam);
    }

}
