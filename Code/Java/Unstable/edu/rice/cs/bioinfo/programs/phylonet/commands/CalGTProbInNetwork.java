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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.*;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
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
    private List<NetworkNonEmpty> _geneTrees;
    private double _bootstrap = 100;

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
        return 10;
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        _speciesNetwork = this.assertAndGetNetwork(0);
        noError = noError && _speciesNetwork != null;

        ParameterIdentList geneTreeParam = this.assertParameterIdentList(1);
        noError = noError && geneTreeParam != null;
        _geneTrees = new LinkedList<NetworkNonEmpty>();
        for(String ident : geneTreeParam.Elements)
        {
            noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
            if(noError)
            {
                _geneTrees.add(this.sourceIdentToNetwork.get(ident));
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

        noError = noError && checkForUnknownSwitches("m", "a", "b");
        checkAndSetOutFile(aParam, mParam, bParam);

        return  noError;
    }

    @Override
    protected String produceResult() {
        StringBuffer result = new StringBuffer();
        
        List<Tree> nbGeneTrees = new ArrayList<Tree>();
        List<Integer> nbCounter = new ArrayList<Integer>();
        for(NetworkNonEmpty geneTree : _geneTrees){
            String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
            NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
            STITree<Double> newtr = new STITree<Double>(true);
            if(_bootstrap<100){
                if(Trees.handleBootStrapInTree(newtr, _bootstrap)==-1){
                    throw new IllegalArgumentException("Input gene tree " + newtr + " have nodes that don't have bootstrap value");
                }

            }
            try
            {
                nr.readTree(newtr);
            }
            catch(Exception e)
            {
                errorDetected.execute(e.getMessage(),
                        this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
            }
            boolean found = false;
            int index = 0;
            for(Tree tr: nbGeneTrees){
                if(Trees.haveSameRootedTopology(tr, newtr)){
                    found = true;
                    break;
                }
                index++;
            }
            if(found){
                nbCounter.set(index, nbCounter.get(index)+1);
            }
            else{
                nbGeneTrees.add(newtr);
                nbCounter.add(1);
            }
        }

        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        Network speciesNetwork = transformer.makeNetwork(_speciesNetwork);

        List<Tree> bGeneTrees = new ArrayList<Tree>();
        List<List<Integer>> nbTree2bTrees = new ArrayList<List<Integer>>();
        for(Tree nbgt: nbGeneTrees){
            List<Integer> bTrees = new ArrayList<Integer>();
            for(Tree bgt: Trees.getAllBinaryResolution(nbgt)){
                int index = 0;
                for(Tree exBgt: bGeneTrees){
                    if(Trees.haveSameRootedTopology(bgt,exBgt)){
                        break;
                    }
                    index++;
                }
                if(index==bGeneTrees.size()){
                    bGeneTrees.add(bgt);
                }
                bTrees.add(index);
            }
            nbTree2bTrees.add(bTrees);
        }

        List<Double> probList;
        if(_multree){
            GeneTreeProbability gtp = new GeneTreeProbability();
            probList = gtp.calculateGTDistribution(speciesNetwork, bGeneTrees, _taxonMap, false);
        }
        else{
            GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
            probList = gtp.calculateGTDistribution(speciesNetwork, bGeneTrees, _taxonMap);
        }
        Iterator<Integer> nbCounterIt = nbCounter.iterator();
        Iterator<List<Integer>> bGTIDs = nbTree2bTrees.iterator();
        double total = 0;
        for(Tree nbgt: nbGeneTrees){
            for(TNode node: nbgt.getNodes()){
                node.setParentDistance(TNode.NO_DISTANCE);
            }
            double maxProb = 0;
            for(int id: bGTIDs.next()){
                maxProb = Math.max(maxProb, probList.get(id));
            }
            int count = nbCounterIt.next();
            total += Math.log(maxProb)*count;
            result.append("\n[x" + count + "] " + nbgt.toString() + " : " + maxProb);
        }
        result.append("\n" + "Total log probability: " + total);

        return result.toString();

    }
}
