/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterQuote;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCURInference_ILP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 12/2/11
 * Time: 1:16 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("inferst_mdc_ilp")
public class InferST_MDC_ILP extends InferSTBase
{

     private Map<String,String> _taxonMap;

    private String _gurobiPath;

    InferST_MDC_ILP(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                    Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
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
        this.richNewickGenerated(tree);

        result.append("\n" + tree +" "+sol._totalCoals+" extra lineages in total");

        return  result.toString();


    }
}
