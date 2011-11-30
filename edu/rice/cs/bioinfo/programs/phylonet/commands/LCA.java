package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.SingleLinePrinter;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.lca.SchieberVishkinLCA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/11
 * Time: 3:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class LCA extends CommandBaseFileOut {

    private NetworkNonEmpty _tree;

    private ArrayList<Set<String>> _setFamilyList = new ArrayList<Set<String>>();

    LCA(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
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

        _tree = this.assertAndGetNetwork(0);
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
