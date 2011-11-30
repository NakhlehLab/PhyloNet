package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.CoalescentCounter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/2/11
 * Time: 2:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class CountCoal extends CommandBaseFileOut {

    private NetworkNonEmpty _tree1, _tree2;

    CountCoal(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 3;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean  noError = true;

        _tree1 = this.assertAndGetNetwork(0);
        noError = noError && _tree1 != null;

        _tree2 = this.assertAndGetNetwork(1);
        noError = noError && _tree2 != null;

        if(this.params.size() == 3)
        {
            noError = noError && checkOutFileContext(2);
        }

        return noError;

    }

    @Override
    protected String produceResult() {

        String eNewickST = NetworkTransformer.toENewick(_tree1);
        String eNewickGT = NetworkTransformer.toENewick(_tree2);

        try
        {
            STITree<Object> st = new STITree<Object>(eNewickST);
		    STITree<Object> gt = new STITree<Object>(eNewickGT);

		    CoalescentCounter cc = new CoalescentCounter();
		    int count = cc.countCoalescents(st, gt);

            return "\n" + count;
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }



    }

}
