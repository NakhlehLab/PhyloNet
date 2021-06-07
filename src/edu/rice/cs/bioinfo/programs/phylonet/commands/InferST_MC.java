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
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MajorityConsensusInference;
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
 * Date: 11/15/11
 * Time: 1:46 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("Infer_ST_MC")
public class InferST_MC extends InferSTBase
{
    private boolean _treesRooted = true;
    private Map<String,String> _taxonMap;
    private int _percentage = 0;

    public InferST_MC(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                      Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 7;  //To change body of implemented methods use File | Settings | File Templates.
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

        ParamExtractor pParam = new ParamExtractor("p", this.params, this.errorDetected);
        if(pParam.ContainsSwitch){
            if(pParam.PostSwitchParam != null)
            {
                try
                {
                    _percentage = Integer.parseInt(pParam.PostSwitchValue);
                    if(_percentage <50){
                        throw new NumberFormatException();
                    }
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized accepted percentage " + pParam.PostSwitchValue, pParam.PostSwitchParam.getLine(), pParam.PostSwitchParam.getColumn());
                }
            }
            else
            {
                errorDetected.execute("Expected value from 50 to 100 after switch -p.", pParam.SwitchParam.getLine(), pParam.SwitchParam.getColumn());
            }
        }

        noError = noError && checkForUnknownSwitches("u", "a", "p");

        this.checkAndSetOutFile(pParam, uParam);

       return  noError;
    }

    @Override
    protected String produceResult() {

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


        MajorityConsensusInference inference = new MajorityConsensusInference();

        Tree inferredTree;
        if(_taxonMap == null){
			inferredTree = inference.inferSpeciesTree(tuples, _treesRooted, _percentage);

		}
		else{
			inferredTree = inference.inferSpeciesTree(tuples,_treesRooted,_taxonMap, _percentage);

		}

        String tree = inferredTree.toStringWD();
        this.richNewickGenerated(tree);
        result.append("\n" + tree);

        return result.toString();
    }
}
