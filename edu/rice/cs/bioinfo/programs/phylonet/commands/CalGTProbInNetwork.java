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

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import sun.nio.cs.ext.MacHebrew;

import java.io.StringReader;
import java.io.StringWriter;
import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 3/13/12
 * Time: 10:56 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("calgtprob")
public class CalGTProbInNetwork extends CommandBaseFileOut{
    private HashMap _taxonMap = null;
    private boolean  _multree = false;
    private NetworkNonEmpty _speciesNetwork;
    private List<List<NetworkNonEmpty>> _geneTrees;
    private int _maxRounds = 100;
    private int _maxTryPerBranch = 100;
    private double _maxBranchLength = 6;
    private double _improvementThreshold = 0.001;
    private double _Brent1 = 0.01;
    private double _Brent2 = 0.001;
    private int _parallel = 1;
    private double _bootstrap = 100;
    private boolean _inferBL = false;
    private int _numRuns = 5;
    private boolean _usingBL = false;
    private boolean _oneGTPerLocus = true;

    public CalGTProbInNetwork(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                              Map<String,NetworkNonEmpty>  sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected,
                              RichNewickReader<Networks> rnReader){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
        return 21;
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

        ParamExtractor mParam = new ParamExtractor("m", this.params, this.errorDetected);
        if(mParam.ContainsSwitch)
        {
            String method = mParam.PostSwitchValue.toLowerCase();
            boolean methodCorrect = method.equals("ac") || method.equals("mul");
            if(!methodCorrect){
                this.errorDetected.execute("-m must be ac or mul", mParam.PostSwitchParam.getLine(), mParam.PostSwitchParam.getColumn());
            }
            else{
                if(method.equals("mul")){
                    _multree = true;
                }
            }
            noError = noError && methodCorrect;
        }

        ParamExtractor aParam = new ParamExtractor("a", this.params, this.errorDetected);
        if(aParam.ContainsSwitch){
            if(_multree){
                ParamExtractorAllelMap aaParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);
                noError = noError && aaParam.IsValidMap;
                if(aaParam.IsValidMap){
                    _taxonMap = aaParam.ValueMap;
                }
            }
            else{
                ParamExtractorAllelListMap aaParam = new ParamExtractorAllelListMap("a", this.params, this.errorDetected);
                noError = noError && aaParam.IsValidMap;
                if(aaParam.IsValidMap){
                    _taxonMap = aaParam.ValueMap;
                }
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

        ParamExtractor rParam = new ParamExtractor("r", this.params, this.errorDetected);
        if(rParam.ContainsSwitch)
        {
            if(rParam.PostSwitchParam != null)
            {
                try
                {
                    _maxRounds = Integer.parseInt(rParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized value after -r " + rParam.PostSwitchValue, rParam.PostSwitchParam.getLine(), rParam.PostSwitchParam.getColumn());
                }
            }
            else
            {
                errorDetected.execute("Expected value after switch -r.", rParam.SwitchParam.getLine(), rParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor tParam = new ParamExtractor("t", this.params, this.errorDetected);
        if(tParam.ContainsSwitch)
        {
            if(tParam.PostSwitchParam != null)
            {
                try
                {
                    _maxTryPerBranch = Integer.parseInt(tParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized value after -t " + tParam.PostSwitchValue, tParam.PostSwitchParam.getLine(), tParam.PostSwitchParam.getColumn());
                }
            }
            else
            {
                errorDetected.execute("Expected value after switch -t.", tParam.SwitchParam.getLine(), tParam.SwitchParam.getColumn());
            }
        }


        ParamExtractor lParam = new ParamExtractor("l", this.params, this.errorDetected);
        if(lParam.ContainsSwitch)
        {
            if(lParam.PostSwitchParam != null)
            {
                try
                {
                    _maxBranchLength = Double.parseDouble(lParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized value of maximum branch length " + lParam.PostSwitchValue, lParam.PostSwitchParam.getLine(), lParam.PostSwitchParam.getColumn());
                }
            }
            else
            {
                errorDetected.execute("Expected value after switch -l.", lParam.SwitchParam.getLine(), lParam.SwitchParam.getColumn());
            }
        }


        ParamExtractor iParam = new ParamExtractor("i", this.params, this.errorDetected);
        if(iParam.ContainsSwitch)
        {
            if(iParam.PostSwitchParam != null)
            {
                try
                {
                    _improvementThreshold = Double.parseDouble(iParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized value after -i " + iParam.PostSwitchValue, iParam.PostSwitchParam.getLine(), iParam.PostSwitchParam.getColumn());
                }
            }
            else
            {
                errorDetected.execute("Expected value after switch -i.", iParam.SwitchParam.getLine(), iParam.SwitchParam.getColumn());
            }
        }


        ParamExtractor pParam = new ParamExtractor("p", this.params, this.errorDetected);
        if(pParam.ContainsSwitch)
        {
            if(pParam.PostSwitchParam != null)
            {
                try
                {
                    if(!(pParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    ParameterIdentList brent = (ParameterIdentList)pParam.PostSwitchParam;
                    int index = 0;
                    for(String value: brent.Elements){
                        if(index==0){
                            _Brent1 = Double.parseDouble(value.trim());
                        }else if(index==1){
                            _Brent2 = Double.parseDouble(value.trim());
                        }
                        else{
                            throw new RuntimeException();
                        }
                        index++;
                    }
                    if(_Brent1==0 || _Brent2==0){
                        throw new RuntimeException();
                    }
                    //ParameterIdentList brent = this.assertParameterIdentList(pParam.PostSwitchValue);
                    /*
                    String brentPrecision = pParam.PostSwitchValue;
                    if(!brentPrecision.contains("~")){
                        throw new NumberFormatException();
                    }
                    String[] splitValues = brentPrecision.split("~");

                    _Brent1 = Double.parseDouble(splitValues[0].trim());
                    _Brent2 = Double.parseDouble(splitValues[1].trim());
                    assert(_Brent1 > _Brent2);
                    */
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized value after switch -p.", pParam.PostSwitchParam.getLine(), pParam.PostSwitchParam.getColumn());
                }

            }
            else
            {
                errorDetected.execute("Expected value after switch -p.", pParam.SwitchParam.getLine(), pParam.SwitchParam.getColumn());
            }
        }


        ParamExtractor plParam = new ParamExtractor("pl", this.params, this.errorDetected);
        if(plParam.ContainsSwitch){
            if(plParam.PostSwitchParam != null)
            {
                try
                {
                    _parallel = Integer.parseInt(plParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized maximum number of networks for search " + plParam.PostSwitchValue, plParam.PostSwitchParam.getLine(), plParam.PostSwitchParam.getColumn());
                }
            }
            else
            {
                errorDetected.execute("Expected value after switch -pl.", plParam.SwitchParam.getLine(), plParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor xParam = new ParamExtractor("x", this.params, this.errorDetected);
        if(xParam.ContainsSwitch)
        {
            if(xParam.PostSwitchParam != null)
            {
                try
                {
                    _numRuns = Integer.parseInt(xParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unrecognized value after -x " + xParam.PostSwitchValue, xParam.PostSwitchParam.getLine(), xParam.PostSwitchParam.getColumn());
                }
            }
            else
            {
                errorDetected.execute("Expected value after switch -x.", xParam.SwitchParam.getLine(), xParam.SwitchParam.getColumn());
            }
        }


        ParamExtractor blParam = new ParamExtractor("bl", this.params, this.errorDetected);
        if(blParam.ContainsSwitch)
        {
            _usingBL = true;
        }


        ParamExtractor oParam = new ParamExtractor("o", this.params, this.errorDetected);
        if(oParam.ContainsSwitch)
        {
            _inferBL = true;
        }


        noError = noError && checkForUnknownSwitches("m", "a", "b", "r", "t", "l", "i", "p", "pl", "x", "bl", "o");
        checkAndSetOutFile(aParam, mParam, bParam, rParam, tParam, lParam, iParam, pParam, plParam, xParam, blParam, oParam);

        return  noError;
    }

    @Override
    protected String produceResult() {
        StringBuffer result = new StringBuffer();

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

        if(!_inferBL){
            _numRuns = 1;
        }
        double optimalScore = Double.NEGATIVE_INFINITY;
        Network optimalNetwork = null;
        if (_usingBL) {
            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            optimalNetwork = transformer.makeNetwork(_speciesNetwork);
            NetworkLikelihoodFromGTTBL scoring;
            if(_oneGTPerLocus){
                scoring = new NetworkLikelihoodFromGTTBL_SingleTreePerLocus();
            }
            else{
                scoring = new NetworkLikelihoodFromGTTBL_MultiTreesPerLocus();
            }
            //ScoreGivenNetworkUsingBLProbabilisticallyParallel scoring = new ScoreGivenNetworkUsingBLProbabilisticallyParallel();
            scoring.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _parallel);
            //optimalScore = scoring.calMLOfNetwork(optimalNetwork, gts, _taxonMap, _inferBL);
            optimalScore = scoring.computeLikelihood(optimalNetwork, gts, _taxonMap, _inferBL);
        }
        else{
            if (_multree) {
                NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
                optimalNetwork = transformer.makeNetwork(_speciesNetwork);

                List<Tree> bGeneTrees = new ArrayList<Tree>();
                List<MutableTuple<List<Integer>,Double>> nbTree2bTrees = new ArrayList<MutableTuple<List<Integer>,Double>>();
                for (List<MutableTuple<Tree, Double>> list : gts) {
                    for (MutableTuple<Tree, Double> nbgt : list) {
                        List<Integer> bTrees = new ArrayList<Integer>();
                        for (Tree bgt : Trees.getAllBinaryResolution(nbgt.Item1)) {
                            int index = 0;
                            for (Tree exBgt : bGeneTrees) {
                                if (Trees.haveSameRootedTopology(bgt, exBgt)) {
                                    break;
                                }
                                index++;
                            }
                            if (index == bGeneTrees.size()) {
                                try {
                                    NewickReader nr = new NewickReader(new StringReader(bgt.toString()));
                                    bGeneTrees.add(nr.readTree());
                                } catch (Exception e) {
                                }

                            }
                            bTrees.add(index);
                        }
                        nbTree2bTrees.add(new MutableTuple<List<Integer>, Double>(bTrees, nbgt.Item2));
                    }
                }
                GeneTreeProbability gtp = new GeneTreeProbability();
                List<Double> scores = gtp.calculateGTDistribution(optimalNetwork, bGeneTrees, _taxonMap, false);

                optimalScore = 0;
                for(MutableTuple<List<Integer>,Double> tuple: nbTree2bTrees) {
                    double subTotal = 0;
                    for (int id : tuple.Item1) {
                        subTotal += scores.get(id);
                    }
                    optimalScore += Math.log(subTotal) * tuple.Item2;
                }

            } else {

                for (int i = 0; i < _numRuns; i++) {
                    NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
                    Network speciesNetwork = transformer.makeNetwork(_speciesNetwork);
                    NetworkLikelihoodFromGTT scoring;
                    if(_oneGTPerLocus){
                        scoring = new NetworkLikelihoodFromGTT_SingleTreePerLocus();
                    }
                    else{
                        scoring = new NetworkLikelihoodFromGTT_MultiTreesPerLocus();
                    }
                    scoring.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _parallel);
                    double score = scoring.computeLikelihood(speciesNetwork, gts, _taxonMap, _inferBL);

                    if (optimalScore < score) {
                        optimalScore = score;
                        optimalNetwork = speciesNetwork;
                    }
                }
            }
        }

        result.append("\nSpecies Network:");
        for(Object node : optimalNetwork.bfs())
        {
            NetNode netNode = (NetNode)node;
            if(!netNode.isLeaf())
            {
                netNode.setName(NetNode.NO_NAME);
            }
        }

        result.append("\n" + optimalNetwork.toString());
        result.append("\n" + "Total log probability: " + optimalScore);
        //result.append(totalScore + " ");

        return result.toString();

    }
}
