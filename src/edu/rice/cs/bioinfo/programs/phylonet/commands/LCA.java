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
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.DAGFactory;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.SingleLinePrinter;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.GraphBuilder;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.lca.SchieberVishkinLCA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/11
 * Time: 3:16 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("LCA")
public class LCA extends CommandBaseFileOut {

    private NetworkNonEmpty _tree;

    private ArrayList<Set<String>> _setFamilyList = new ArrayList<Set<String>>();

    public LCA(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    protected int getMinNumParams()
    {
        return 2;
    }

    protected int getMaxNumParams()
    {
        return 3;
    }

    public boolean checkParamsForCommand() {

        boolean noError = true;

        _tree = this.assertAndGetTree(0);
        noError = noError && _tree != null;


        final Parameter setFamilyParam = params.get(1);

        noError = setFamilyParam.execute(new ParameterAlgo<Boolean, Boolean, RuntimeException>() {
            public Boolean forIdentifier(ParameterIdent parameter, Boolean o) throws RuntimeException {
                return unExpectedCase();
            }

            public Boolean forIdentList(ParameterIdentList parameterIdentList, Boolean aBoolean) throws RuntimeException {
                return unExpectedCase();
            }

            public Boolean forQuote(ParameterQuote parameter, Boolean o) throws RuntimeException {
                return unExpectedCase();
            }

            public Boolean forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Boolean noError) throws RuntimeException {

                for(Iterable<String> set : parameterTaxonSetList.TaxonSetList)
                {
                    HashSet<String> taxa = new HashSet<String>();
                    for(String element : set)
                    {
                        if(!taxa.contains(element))
                        {
                            taxa.add(element);
                        }
                        else
                        {
                            errorDetected.execute(
                                    String.format("Duplicate taxon '%s'.", element), setFamilyParam.getLine(), setFamilyParam.getColumn());
                            noError = false;
                        }
                    }
                    _setFamilyList.add(taxa);
                }

                return noError;
            }

            public Boolean forIdentSet(ParameterIdentSet parameterIdentSet, Boolean o) throws RuntimeException {
                return unExpectedCase();
            }

            public Boolean forTaxaMap(ParameterTaxaMap parameterTaxaMap, Boolean aBoolean) throws RuntimeException {
                return unExpectedCase();
            }

            private Boolean unExpectedCase()
            {
                errorDetected.execute("Expected a set list. (E.g. '({A,B},{C,D})')", setFamilyParam.getLine(), setFamilyParam.getColumn());
                return false;
            }
        }, noError);

        if(params.size() == 3)
        {
            noError = this.checkOutFileContext(2);
        }

        /**
         * Check to see all the taxon names in the set family appear in the tree.
         */
        if(noError)
        {
            final HashSet<String> taxonNames = new HashSet<String>();
            DAGFactory.makeDAG(_tree, new GraphBuilder<Object>() {
                public Object createNode(String s) {
                    if (!taxonNames.contains(s))
                        taxonNames.add(s);
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Object createHybridNode(String s, HybridNodeType hybridNodeType, BigInteger bigInteger) {
                    // its a tree, so this should never be called.
                    return null;
                }

                public void createDirectedEdge(Object o, Object o1, BigDecimal bigDecimal, BigDecimal bigDecimal1, BigDecimal bigDecimal2) {
                    //To change body of implemented methods use File | Settings | File Templates.
                }
            });

            for(Set<String> set : _setFamilyList)
            {
                for(String taxonName : set)
                {
                    if(!taxonNames.contains(taxonName))
                    {
                        ParameterIdent treeIdent = (ParameterIdent)this.params.get(0);
                        this.errorDetected.execute(String.format("Taxon name '%s' does not appear in tree '%s'.", taxonName, treeIdent.Content), setFamilyParam.getLine(), setFamilyParam.getColumn());
                        noError = false;
                    }
                }
            }
        }

        return noError;
    }

    private boolean addTaxonToSet(StringBuffer taxonNameUnderConstruction, HashSet<String> set, Proc3<String,Integer,Integer> errorDetected, int lineNum, int i, Parameter setFamilyParam)
    {
        String taxonName = taxonNameUnderConstruction.toString().trim();
        if(set.contains(taxonName))
        {
            errorDetected.execute(String.format("Duplicate taxon '%s'.'", taxonName), lineNum, setFamilyParam.getColumn() + i + 1 );
            return false;
        }
        else
        {
            set.add(taxonName);
            return true;
        }
    }

    public String produceResult()
    {

        StringBuilder result = new StringBuilder();

        MutableTree tree = NetworkTransformer.toTree(_tree);
        Trees.autoLabelNodes(tree);

        NetworkNonEmpty autoLabeledNetwork = (NetworkNonEmpty) TreeTransformer.toNetwork(tree);
        String richNewickString = new SingleLinePrinter().toString(autoLabeledNetwork);

        result.append("\n" + richNewickString);
        this.richNewickGenerated(richNewickString);


        SchieberVishkinLCA solver = new SchieberVishkinLCA(tree);



        for(Set<String> taxaSet : _setFamilyList)
        {
            TNode lca = null;

            for(String taxon : taxaSet)
            {
                if(lca == null)
                {
                    lca = tree.getNode(taxon);

                    if(lca == null)
                    {
                        throw new RuntimeException("Taxon not a member of tree.");
                    }
                    continue;
                }

                TNode node = tree.getNode(taxon);

                if(node == null)
                {
                    throw new RuntimeException("Taxon not a member of tree.");
                }

                lca = solver.getLCA(lca, node);
            }

            result.append("\n" + lca.getName());
        }

        return result.toString();

    }

}
