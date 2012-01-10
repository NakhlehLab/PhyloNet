package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterQuote;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.DeepCoalescencesCounter;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCURInference_ILP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeClusterWD;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.awt.geom.Path2D;
import java.io.*;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 12/2/11
 * Time: 1:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferST_MDC_ILP extends InferSTBase
{

     private Map<String,String> _taxonMap;

    private String _gurobiPath;

    InferST_MDC_ILP(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 5;  //To change body of implemented methods use File | Settings | File Templates.
    }


    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = super.checkParamsForCommand();

        ParameterQuote gurobiPath = this.assertQuotedParameter(1);
        noError = noError && gurobiPath != null;

        if(gurobiPath != null)
        {
            if(!new File(gurobiPath.UnquotedText).exists())
            {
                this.errorDetected.execute("Invalid gurobi gurobi path: '" + _gurobiPath + "'.", gurobiPath.getLine(), gurobiPath.getColumn());
                noError = false;
            }
            else
            {
                _gurobiPath = gurobiPath.UnquotedText;
            }
        }


        TaxonMapResult result = assignTaxonMap();
        noError = noError && result.NoError;
        _taxonMap = result.TaxonMap;

         noError = noError && checkForUnknownSwitches("a");

        this.checkAndSetOutFile(result.Extractor);

       return  noError;
    }

    @Override
    protected String produceResult() {

        StringBuffer result = new StringBuffer();

        MDCURInference_ILP inference = new MDCURInference_ILP();
       Solution sol;

        List<Tree> trees = GetGeneTreesAsTreeList();

        if(_taxonMap == null){
			sol = inference.inferSpeciesTree(_gurobiPath,trees);
		}
		else {
			sol = inference.inferSpeciesTree(_gurobiPath,trees, _taxonMap);
		}

        String tree = sol._st.toStringWD();
        this.speciesTreeGenerated(tree);

        result.append("\n" + tree +" "+sol._totalCoals+" extra lineages in total");

        return  result.toString();


    }
}
