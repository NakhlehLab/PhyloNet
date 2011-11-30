package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MajorityConsensusInference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/15/11
 * Time: 1:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferST_MC extends InferSTBase
{
    private boolean _treesRooted = true;
    private Map<String,String> _taxonMap;

    InferST_MC(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 5;  //To change body of implemented methods use File | Settings | File Templates.
    }


    @Override
    protected boolean checkParamsForCommand() {

       boolean noError = super.checkParamsForCommand();

         TaxonMapResult result = assignTaxonMap();
        noError = noError && result.NoError;
        _taxonMap = result.TaxonMap;

        ParamExtractor uParam = new ParamExtractor("u", this.params, this.errorDetected);
        if(uParam.ContainsSwitch)
        {
            _treesRooted = false;
        }

        this.checkAndSetOutFile();

       return  noError;
    }

    @Override
    protected String produceResult() {

        if(_geneTrees == null)
        {
            throw new IllegalStateException();
        }

        StringBuffer result = new StringBuffer();

        List<Tree> trees = GetGeneTreesAsSTIDoubleTreeList();
        MajorityConsensusInference inference = new MajorityConsensusInference();

        Tree inferredTree;
        if(_taxonMap == null){
			if(_treesRooted){
				inferredTree = inference.inferSpeciesTreeRooted(trees);
			}
			else{
				inferredTree = inference.inferSpeciesTreeUnrooted(trees);
			}
		}
		else{
			if(_treesRooted){
				inferredTree = inference.inferSpeciesTreeRooted(trees,_taxonMap);
			}
			else {
				inferredTree = inference.inferSpeciesTreeUnrooted(trees,_taxonMap);
			}
		}

        result.append("\n" + inferredTree.toStringWD());

        return result.toString();
    }
}
