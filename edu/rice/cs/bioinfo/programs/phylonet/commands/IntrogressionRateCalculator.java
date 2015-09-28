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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterTaxonSetList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.introgression.CoalescenceAnalysis;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: dw20
 * Date: 9/28/15
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("CalIntroRate")
public class IntrogressionRateCalculator extends CommandBaseFileOut{


    private NetworkNonEmpty _speciesNetwork;
    private List<List<NetworkNonEmpty>> _geneTrees;
    private boolean _oneGTPerLocus = true;
    private double _bootstrap = 100;
    private boolean _xlg = false;

    public IntrogressionRateCalculator(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                                       Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected,
                                       RichNewickReader<Networks> rnReader){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
        return 4;
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        _speciesNetwork = this.assertAndGetNetwork(0);
        noError = noError && _speciesNetwork != null;

        Parameter geneTreeParam = this.assertParameterIdentListOrSetList(1);
        noError = noError && geneTreeParam != null;
        _geneTrees = new LinkedList<List<NetworkNonEmpty>>();

        if(geneTreeParam instanceof ParameterIdentList) {
            ParameterIdentList geneTreeList = (ParameterIdentList)geneTreeParam;
            for (String ident : geneTreeList.Elements) {
                noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
                if (noError) {
                    _geneTrees.add(Arrays.asList(this.sourceIdentToNetwork.get(ident)));
                }
            }
        }
        else{
            _oneGTPerLocus = false;
            ParameterTaxonSetList geneTreeSetList = (ParameterTaxonSetList)geneTreeParam;
            for(Iterable<String> gtSet : geneTreeSetList.TaxonSetList)
            {
                List<NetworkNonEmpty> geneTreesForOneLocus = new ArrayList<NetworkNonEmpty>();
                for (String ident : gtSet) {
                    noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
                    if (noError) {
                        geneTreesForOneLocus.add(this.sourceIdentToNetwork.get(ident));
                    }
                }
                _geneTrees.add(geneTreesForOneLocus);
            }
        }

        ParamExtractor bParam = new ParamExtractor("b", this.params, this.errorDetected);
        if(bParam.ContainsSwitch)
        {
            if(bParam.PostSwitchParam != null)
            {
                try
                {
                    _bootstrap = Double.parseDouble(bParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized bootstrap value " + bParam.PostSwitchValue, bParam.PostSwitchParam.getLine(), bParam.PostSwitchParam.getColumn());
                }
            }
            else
            {
                errorDetected.execute("Expected value after switch -b.", bParam.SwitchParam.getLine(), bParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor lParam = new ParamExtractor("l", this.params, this.errorDetected);
        if(lParam.ContainsSwitch)
        {
            _xlg = true;
        }


        noError = noError && checkForUnknownSwitches("b") && checkForUnknownSwitches("l");
        checkAndSetOutFile(bParam, lParam);

        return  noError;
    }

    @Override
    protected String produceResult() {

        // gene trees
        List<List<MutableTuple<Tree,Double>>> gts = new ArrayList<List<MutableTuple<Tree,Double>>>();
        for(List<NetworkNonEmpty> geneTrees : _geneTrees) {
            List<MutableTuple<Tree, Double>> gtsForOneLocus = new ArrayList<MutableTuple<Tree, Double>>();
            for (NetworkNonEmpty geneTree : geneTrees) {
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

                try {
                    nr.readTree(newtr);
                } catch (Exception e) {
                    errorDetected.execute(e.getMessage(),
                            this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
                }
                Trees.removeBinaryNodes(newtr);
                if (_bootstrap < 100) {
                    if (Trees.handleBootStrapInTree(newtr, _bootstrap) == -1) {
                        throw new IllegalArgumentException("Input gene tree " + newtr + " have nodes that don't have bootstrap value");
                    }

                }
                gtsForOneLocus.add(new MutableTuple<Tree, Double>(newtr, weight));

            }
            gts.add(gtsForOneLocus);
        }
        // compute
        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        Network net = transformer.makeNetwork(_speciesNetwork);
        CoalescenceAnalysis analysis = new CoalescenceAnalysis(net, gts);
        return _xlg ? analysis.computeXLG() : analysis.computeIntrogressionRate();
    }
}
