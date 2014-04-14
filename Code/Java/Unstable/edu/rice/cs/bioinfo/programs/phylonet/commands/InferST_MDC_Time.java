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
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCWithTime;
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
 * Date: 11/9/11
 * Time: 4:54 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("Infer_ST_MDC_Time")
public class InferST_MDC_Time extends InferSTBase
{

    private Map<String,String> _taxonMap;
	private double _bootstrap = 100;


    public InferST_MDC_Time(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                            Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 6;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = super.checkParamsForCommand();


        ThresholdResult tr = this.assignThreshold(100);
        _bootstrap = tr.Threshold;
        noError = noError && tr.NoError;

         TaxonMapResult result = assignTaxonMap();
        noError = noError && result.NoError;
        _taxonMap = result.TaxonMap;

         noError = noError && checkForUnknownSwitches("a","b");

        this.checkAndSetOutFile(tr.Extractor, result.Extractor);

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


        MDCWithTime inference = new MDCWithTime();

		Solution s;
		if(_taxonMap == null){
			s = inference.inferSpeciesTree(tuples, _bootstrap);
		}
		else{
			s = inference.inferSpeciesTree(tuples,_taxonMap, _bootstrap);
		}

        String tree = s._st.toStringWD();
        this.richNewickGenerated(tree);
        result.append("\n" + s._st.toStringWD()+" "+s._totalCoals+" extra lineages in total");

        return result.toString();
    }

}
