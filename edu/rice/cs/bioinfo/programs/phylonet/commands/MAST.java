package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.ContainsHybridNode;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.SingleLinePrinter;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.mast.SteelWarnowMAST;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.ExMultipleMasts;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/29/11
 * Time: 2:00 PM
 * To change this template use File | Settings | File Templates.
 */
class MAST extends CommandBaseFileOut
{
    private boolean _computeAll = false;

   // Parameter _potentialTreeSetParam = null;

  //  ParameterIdentSet _treeSetParameter = null;

    Set<Tree> _treeSet = new HashSet<Tree>();

    boolean _allUnrooted = true;

    MAST(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    protected int getMinNumParams()
    {
        return 1;
    }

     protected int getMaxNumParams()
    {
        return 3;
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

                    if(rootedness != null)
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

        final String allSwitchText = "-a";
        ParamExtractor allSwitch = new ParamExtractor("a", this.params, this.errorDetected);
        _computeAll = allSwitch.ContainsSwitch;

        if(params.size() == 3)
        {
            noError = noError && this.checkOutFileContext(2);

            if(!allSwitch.ContainsSwitch)
            {
                Parameter param0 = this.params.get(0);
              errorDetected.execute(
                            String.format("Expected '%s' switch for command '%s'.", allSwitchText, this.getDefiningSyntaxCommand()),
                            param0.getLine(), param0.getColumn());
                    noError = false;
            }
        }
        else if(!allSwitch.ContainsSwitch && params.size() == 2)
        {
            noError = noError && this.checkOutFileContext(1);
        }

         noError = noError && checkForUnknownSwitches("a");

        if(!noError)
        {
            _treeSet = null;
        }

        return noError;

    }


    @Override
    protected String produceResult() {

        StringBuilder result = new StringBuilder();

        if(_computeAll)
        {
            ExMultipleMasts emm = new ExMultipleMasts();

            Iterator<Tree> trees = _treeSet.iterator();
            Set<Tree> mastSet = emm.computeMultipleMasts(trees.next(), trees.next());

            result.append("\nNumber of MASTs: " + mastSet.size());

            for (Tree mast : mastSet) {
                NetworkNonEmpty netMast = (NetworkNonEmpty) TreeTransformer.toNetwork(mast);
                String netMastString = new SingleLinePrinter().toString(netMast);
                this.treeGenerated(netMastString);
				result.append("\n" + netMastString);
			}
        }
        else
        {
            SteelWarnowMAST calculator = new SteelWarnowMAST();
            Tree outTree;

            if(_allUnrooted)
            {
                outTree = calculator.computeUMAST(_treeSet);
            }
            else
            {
                outTree = calculator.computeRMAST(_treeSet);
            }

            NetworkNonEmpty netMast = (NetworkNonEmpty) TreeTransformer.toNetwork(outTree);
            String netMastString = new SingleLinePrinter().toString(netMast);
            this.treeGenerated(netMastString);
            result.append("\n" + netMastString);
        }

        return result.toString();

    }
}
