package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test.multi;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test.HCGModelBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variationalMulti.VariationalModelMulti;


public class TestMultiModel {
    public static void test1() {
        ModelTree model = HCGModelBuilder.getHCGModel();
        VariationalModelMulti variationalPosterior = new VariationalModelMulti(model);
    }

    public static void main(String[] args) {
        test1();
    }
}
