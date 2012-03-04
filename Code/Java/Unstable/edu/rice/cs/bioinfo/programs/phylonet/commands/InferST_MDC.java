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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/9/11
 * Time: 2:42 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferST_MDC extends InferSTBase
{

    private Map<String,String> _taxonMap;

    private boolean _explore = false;
    private double _proportion = 0;
	private boolean _exhaust = false;
	private double _bootstrap = 1;
	private double _time = -1;
	private boolean _unresolved = false;

    InferST_MDC(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 14;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = super.checkParamsForCommand();

        ParamExtractor eParam = new ParamExtractor("e", this.params, this.errorDetected);
        if(eParam.ContainsSwitch)
        {
            if(eParam.PostSwitchParam != null)
            {
                try
                {
                    _proportion = Double.parseDouble(eParam.PostSwitchValue);
                    _explore = true;
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unknown proportion value '" + eParam.PostSwitchValue + "'.", eParam.PostSwitchParam.getLine(), eParam.PostSwitchParam.getColumn());
                    noError = false;
                }
            }
            else
            {
                errorDetected.execute("Expected value after proportion switch.", eParam.SwitchParam.getLine(), eParam.SwitchParam.getColumn());
                noError = false;
            }
        }

        ParamExtractor tParam = new ParamExtractor("t", this.params, this.errorDetected);
        if(tParam.ContainsSwitch)
        {
            if(tParam.PostSwitchParam != null)
            {
                try
                {
                    _time = Double.parseDouble(tParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unknown time value '" + tParam.PostSwitchValue + "'.", tParam.PostSwitchParam.getLine(), tParam.PostSwitchParam.getColumn());
                    noError = false;
                }
            }
            else
            {
                errorDetected.execute("Expected value after time switch.", tParam.SwitchParam.getLine(), tParam.SwitchParam.getColumn());
                noError = false;
            }
        }

        ParamExtractor xParam = new ParamExtractor("x", this.params, this.errorDetected);
        if(xParam.ContainsSwitch)
        {
            _exhaust = true;
        }

        ParamExtractor urParam = new ParamExtractor("ur", this.params, this.errorDetected);
        if(urParam.ContainsSwitch)
        {
            _unresolved = true;
        }



        ThresholdResult tr = this.assignThreshold(1);
        _bootstrap = tr.Threshold;
        noError = noError && tr.NoError;


        TaxonMapResult result = assignTaxonMap();
        noError = noError && result.NoError;
        _taxonMap = result.TaxonMap;

         noError = noError && checkForUnknownSwitches("e", "t", "x", "ur", "b", "a");

        this.checkAndSetOutFile(eParam, tr.Extractor, result.Extractor, tParam);


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

        List<Tree> trees = GetGeneTreesAsSTIDoubleTreeList();

        MDCInference_DP inference = new MDCInference_DP();

		List<Solution> solutions;
		if(_taxonMap == null){
			solutions = inference.inferSpeciesTree(trees,_explore,_proportion,_exhaust,_bootstrap,_unresolved,_time);
		}
		else {
			solutions = inference.inferSpeciesTree(trees,_taxonMap,_explore,_proportion,_exhaust,_bootstrap,_unresolved,_time);
		}

        for(Solution s: solutions){
                String tree = s._st.toStringWD();
                this.treeGenerated(tree);
				result.append("\n" + s._st.toStringWD()+" "+s._totalCoals+" extra lineages in total");
			}

        return result.toString();

    }
}

