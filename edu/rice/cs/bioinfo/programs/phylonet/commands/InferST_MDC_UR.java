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
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCURInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/10/11
 * Time: 1:20 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("Infer_ST_MDC_UR")
public class InferST_MDC_UR extends InferSTBase
{
    private Map<String,String> _taxonMap;
    private boolean _explore = false;
    private double _proportion = 0;
	private boolean _exhaust = false;
	private double _bootstrap = 100;
	private double _time = -1;
	private boolean _unresolved = false;

    public InferST_MDC_UR(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                          Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 12;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand()
    {
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



        ThresholdResult tr = this.assignThreshold(100);
        _bootstrap = tr.Threshold;
        noError = noError && tr.NoError;

        TaxonMapResult result = assignTaxonMap();
        noError = noError && result.NoError;
        _taxonMap = result.TaxonMap;

         noError = noError && checkForUnknownSwitches("e", "t", "x", "ur", "a", "b");

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

        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();

        for(NetworkNonEmpty geneTree : _geneTrees)
        {
            double weight = geneTree.TreeProbability.execute(new TreeProbabilityAlgo<Double, RuntimeException>() {
                @Override
                public Double forEmpty(TreeProbabilityEmpty empty) {
                    return 1.0;
                }

                @Override
                public Double forNonEmpty(TreeProbabilityNonEmpty nonEmpty) {
                    return Double.parseDouble(nonEmpty.ProbString);
                }
            });

            String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
            NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
            STITree<Double> newtr = new STITree<Double>(true);
            try
            {
                nr.readTree(newtr);
            }
            catch(Exception e)
            {
                errorDetected.execute(e.getMessage(),
                        this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
            }
            Trees.removeBinaryNodes(newtr);
            if(_bootstrap<100){
                if(Trees.handleBootStrapInTree(newtr, _bootstrap)==-1){
                    throw new IllegalArgumentException("Input gene tree " + newtr + " have nodes that don't have bootstrap value");
                }

            }
            for(TNode node: newtr.getNodes()){
                node.setParentDistance(TNode.NO_DISTANCE);
            }

            String exp = Trees.getLexicographicNewickString(newtr, _taxonMap);
            MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
            if(existingTuple==null){
                existingTuple = new MutableTuple<Tree, Double>(newtr, weight);
                exp2tree.put(exp, existingTuple);
            }
            else{
                existingTuple.Item2 += weight;
            }
        }

        List<MutableTuple<Tree,Double>> tuples = new ArrayList<MutableTuple<Tree, Double>>();
        tuples.addAll(exp2tree.values());



        MDCURInference_DP inference = new MDCURInference_DP();

        List<Solution> solutions;
		if(_taxonMap == null){
			solutions = inference.inferSpeciesTree(tuples,_explore,_proportion,_exhaust,_unresolved,_time);
		}
		else
            solutions = inference.inferSpeciesTree(tuples,_taxonMap,_explore,_proportion,_exhaust,_unresolved,_time);

			for(Solution s: solutions){
                String tree =  s._st.toStringWD();
                this.richNewickGenerated(tree);
				result.append("\n" + tree +" "+s._totalCoals+" extra lineages in total");
			}


        return result.toString();
    }
}
