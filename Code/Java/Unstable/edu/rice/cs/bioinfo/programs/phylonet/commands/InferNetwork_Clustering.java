package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentSet;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.ContainsHybridNode;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.InferNetworkClustering;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.lang.reflect.Array;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/7/16
 * Time: 1:29 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("InferNetwork_Clustering")
public class InferNetwork_Clustering extends CommandBaseFileOut {

    List<Tree> _treeSet = new ArrayList<>();

    public InferNetwork_Clustering(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected,
                RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    protected int getMinNumParams()
    {
        return 1;
    }

    protected int getMaxNumParams()
    {
        return 1;
    }

    public boolean checkParamsForCommand() {

        boolean noError = true;

        ParameterIdentSet treesIdentSet = this.assertParameterIdentSet(0);
        noError = noError && treesIdentSet != null;

        if(treesIdentSet != null)
        {
            // LinkedList<NetworkNonEmpty> networks = this.assertNetworksExist(treesIdentSet);

            Boolean rootedness = null;
            for(String potentialTreeIdent : treesIdentSet.Elements)
            {
                NetworkNonEmpty potentialTree = this.sourceIdentToNetwork.get(potentialTreeIdent);

                if(potentialTree.execute(ContainsHybridNode.Singleton, null))
                {
                    noError = false;
                    this.errorDetected.execute(String.format("Expected '%s' to be a tree but contains a hybrid node.", potentialTreeIdent)
                            ,treesIdentSet.getLine(), treesIdentSet.getColumn());
                }
                else
                {
                    Tree treeForm = NetworkTransformer.toTree(potentialTree);

                    /*if(rootedness != null)
                    {
                        if(rootedness != treeForm.isRooted())
                        {
                            errorDetected.execute(
                                    String.format("All trees for MAST must be rooted or unrooted."),
                                    treesIdentSet.getLine(), treesIdentSet.getColumn());
                            noError = false;
                        }
                    }
                    else
                    {
                        rootedness = treeForm.isRooted();
                        _allUnrooted = !treeForm.isRooted();
                    }*/

                    if(!treeForm.isRooted()) {
                        errorDetected.execute(
                                String.format("All trees must be rooted."),
                                treesIdentSet.getLine(), treesIdentSet.getColumn());
                        noError = false;
                    }

                    if(treeForm.getLeafCount() < 3)
                    {
                        noError = false;
                        this.errorDetected.execute(String.format("All trees to MAST must have at least three leaves. '%s' does not.", potentialTreeIdent),
                                treesIdentSet.getLine(), treesIdentSet.getColumn());
                    }

                    _treeSet.add(treeForm);
                }
            }


        }

        if(!noError)
            _treeSet = null;

        return noError;
    }

    @Override
    protected String produceResult() {
        StringBuilder result = new StringBuilder();

        System.out.println();

        InferNetworkClustering inferNetworkClustering = new InferNetworkClustering();
        List<List<MutableTuple<Tree, Double>>> gts = new ArrayList<>();
        for(Tree tree : _treeSet) {
            gts.add(Arrays.asList(new MutableTuple<Tree, Double>(tree, 1.0)));
        }
        result.append(inferNetworkClustering.inferNetwork(gts));

        return result.toString();
    }
}
