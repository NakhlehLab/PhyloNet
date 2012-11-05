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

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
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
public class CalGTProbInNetwork extends CommandBaseFileOut{
    private HashMap _taxonMap = null;
    private boolean  _multree = false;
    private NetworkNonEmpty _speciesNetwork;
    private List<NetworkNonEmpty> _geneTrees;

    public CalGTProbInNetwork(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                      Map<String,NetworkNonEmpty>  sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
        return 8;
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

        noError = noError && checkForUnknownSwitches("m", "a");
        checkAndSetOutFile(aParam, mParam);

        return  noError;
    }

    @Override
    protected String produceResult() {
        StringBuffer result = new StringBuffer();
        
        List<Tree> geneTrees = new ArrayList<Tree>();
        List<Integer> counter = new ArrayList<Integer>();
        for(NetworkNonEmpty geneTree : _geneTrees){
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
            boolean found = false;
            int index = 0;
            for(Tree tr: geneTrees){
                if(Trees.haveSameRootedTopology(tr, newtr)){
                    found = true;
                    break;
                }
                index++;
            }
            if(found){
                counter.set(index, counter.get(index)+1);
            }
            else{
                geneTrees.add(newtr);
                counter.add(1);
            }
        }

        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        Network speciesNetwork = transformer.makeNetwork(_speciesNetwork);


        Iterator<Double> probList;
        if(_multree){
            GeneTreeProbability gtp = new GeneTreeProbability();
            probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, _taxonMap, false).iterator();
        }
        else{
            GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
            probList = gtp.calculateGTDistribution(speciesNetwork, geneTrees, _taxonMap).iterator();
        }
        Iterator<Integer> counterIt = counter.iterator();
        double total = 0;
        for(Tree gt: geneTrees){
            for(TNode node: gt.getNodes()){
                node.setParentDistance(TNode.NO_DISTANCE);
            }
            double prob = probList.next();
            int count = counterIt.next();
            total += Math.log(prob)*count;
            result.append("\n[x" + count + "] " + gt.toString() + " : " + prob);
        }
        result.append("\n" + "Total log probability: " + total);

        return result.toString();

    }
}
