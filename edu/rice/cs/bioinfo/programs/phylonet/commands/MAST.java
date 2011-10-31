package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
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
    Parameter _outFile = null;

    Parameter _allSwitch = null;

    Parameter _potentialTreeSetParam = null;

    ParameterIdentSet _treeSetParameter = null;

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

        final String allSwitchText = "-a";

        Parameter potentialTreeSetParam = null;
        if(params.size() == 3)
        {
            _allSwitch = params.get(0);
            _potentialTreeSetParam = params.get(1);
            _outFile = params.get(2);
        }
        else if(params.size() == 2)
        {
            noError = params.get(0).execute(new ParameterAlgo<Boolean, Boolean, RuntimeException>() {

                public Boolean forIdentifier(ParameterIdent parameter, Boolean noError) throws RuntimeException {
                   _allSwitch = params.get(0);
                   _potentialTreeSetParam = params.get(1);
                   return noError;
                }

                public Boolean forIdentList(ParameterIdentList parameterIdentList, Boolean aBoolean) throws RuntimeException {
                   return false;
                }

                public Boolean forQuote(ParameterQuote parameter, Boolean noError) throws RuntimeException {
                   return false;
                }

                public Boolean forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Boolean noError) throws RuntimeException {
                    return false;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Boolean forIdentSet(ParameterIdentSet parameterIdentSet, Boolean noError) throws RuntimeException {
                    _treeSetParameter = parameterIdentSet;
                    return noError;
                }

                public Boolean forTaxaMap(ParameterTaxaMap parameterTaxaMap, Boolean aBoolean) throws RuntimeException {
                    return false;
                }

            }, noError);
        }
        else if(params.size() == 1)
        {
            _potentialTreeSetParam = params.get(0);
        }

        if(_allSwitch != null)
        {
            String switchText = _allSwitch.execute(GetSimpleParamValue.Singleton, null);

            if(switchText != null)
            {
                if(!switchText.toLowerCase().equals(allSwitchText))
                {
                    errorDetected.execute(
                            String.format("Expected '%s' switch for command '%s'.", allSwitchText, this.getDefiningSyntaxCommand()),
                            _allSwitch.getLine(), _allSwitch.getColumn());
                    noError = false;
                }
            }
        }

        if(_treeSetParameter == null && _potentialTreeSetParam != null)
        {
            try
            {
                _treeSetParameter = (ParameterIdentSet)_potentialTreeSetParam;
            }
            catch(ClassCastException e)
            {
                errorDetected.execute("Expected set of networks or trees.",
                        _potentialTreeSetParam.getLine(), _potentialTreeSetParam.getColumn());
                noError = false;
            }
        }


        if(noError)
        {
            return checkContext(sourceIdentToNetwork, errorDetected);
        }
        else
        {
            return noError;
        }

    }

    private boolean checkContext(Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected)
    {
        Boolean noError = true;

        if(_outFile != null)
        {
            noError = this.checkOutFileContext(_outFile, errorDetected);
        }

        Boolean rootedness = null;

        for(String networkName : _treeSetParameter.Elements)
        {
            noError = noError && this.assertNetworkExists(networkName, _treeSetParameter.getLine(), _treeSetParameter.getColumn());

            if(noError)
            {
                Network n = sourceIdentToNetwork.get(networkName);
                Tree treeForm = NetworkTransformer.toTree(n);

                if(rootedness != null)
                {
                   if(rootedness != treeForm.isRooted())
                   {
                       errorDetected.execute(
                               String.format("All trees for for MAST must be rooted or unrooted.  Unexpected rootendess for tree '%s'.", networkName),
                               this.getDefiningSyntaxCommand().getLine(), this.getDefiningSyntaxCommand().getColumn());
                       noError = false;
                   }
                }
                else
                {
                    rootedness = treeForm.isRooted();
                    _allUnrooted = !treeForm.isRooted();
                }

                _treeSet.add(treeForm);
            }
        }

         return  noError;


    }

    @Override
    protected String produceResult() {

        StringBuilder result = new StringBuilder();

        if(_allSwitch != null)
        {
            ExMultipleMasts emm = new ExMultipleMasts();

            Iterator<Tree> trees = _treeSet.iterator();
            Set<Tree> mastSet = emm.computeMultipleMasts(trees.next(), trees.next());

            result.append("Number of MASTs: " + mastSet.size());

            for (Tree mast : mastSet) {
                NetworkNonEmpty netMast = (NetworkNonEmpty) TreeTransformer.toNetwork(mast);
				result.append("\n" + new SingleLinePrinter().toString(netMast));
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
            result.append(new SingleLinePrinter().toString(netMast));
        }

        return result.toString();

    }
}
