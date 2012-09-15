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
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.MDCOnNetwork;
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
public class CountXLInNetwork extends CommandBaseFileOut{
    private HashMap<String,String> _taxonMap = null;
    private NetworkNonEmpty _speciesNetwork;
    private List<NetworkNonEmpty> _geneTrees;

    public CountXLInNetwork(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                              Map<String,NetworkNonEmpty>  sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
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

        ParameterIdentList geneTreeParam = this.assertParameterIdentList(1);
        noError = noError && geneTreeParam != null;
        _geneTrees = new LinkedList<NetworkNonEmpty>();

        if(noError)
        {

            for(String ident : geneTreeParam.Elements)
            {
                noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
                if(noError)
                {
                    _geneTrees.add(this.sourceIdentToNetwork.get(ident));
                }
            }

            ParamExtractorAllelMap aParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);
            if(aParam.ContainsSwitch){
                noError = noError && aParam.IsValidMap;
                if(aParam.IsValidMap){
                    _taxonMap = aParam.ValueMap;
                }
            }

             noError = noError && checkForUnknownSwitches("a");
             checkAndSetOutFile(aParam);
        }



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

        MDCOnNetwork mdc = new MDCOnNetwork();
        Iterator<Integer> xlList = mdc.countExtraCoal(speciesNetwork, geneTrees, _taxonMap).iterator();
        Iterator<Integer> counterIt = counter.iterator();
        int total = 0;
        for(Tree gt: geneTrees){
            for(TNode node: gt.getNodes()){
                node.setParentDistance(TNode.NO_DISTANCE);
            }
            int xl = xlList.next();
            int count = counterIt.next();;
            total += xl * count;
            result.append("\n[x" + count + "] " + gt.toString() + ": " + xl);
        }
        result.append("\n" + "Total number of extra lineages: " + total);

        return result.toString();

    }
}
