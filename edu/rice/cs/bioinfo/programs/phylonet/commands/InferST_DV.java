package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.DemocraticVoteInference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/14/11
 * Time: 4:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferST_DV extends InferSTBase
{
    InferST_DV(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand()
    {
       boolean noError = super.checkParamsForCommand();

        this.checkAndSetOutFile();


       return  noError;
    }

    @Override
    protected String produceResult()
    {
        if(_geneTrees == null)
        {
            throw new IllegalStateException();
        }

        StringBuffer result = new StringBuffer();

        List<Tree> trees = GetGeneTreesAsTreeList();
        DemocraticVoteInference inference = new DemocraticVoteInference();
		List<Tree> inferredTrees = inference.inferSpeciesTree(trees);

        result.append("\n" + inferredTrees.size() + " tree(s) has(have) the frequency of "
					+ inference.getFrequency() + "%");
			for (Tree tr : inferredTrees){
                String tree = tr.toString();
                this.SpeciesTreeGenerated(tree);
				result.append("\n" + tree);
			}

        return result.toString();


    }
}
