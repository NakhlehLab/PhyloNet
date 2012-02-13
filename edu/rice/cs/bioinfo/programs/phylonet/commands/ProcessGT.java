package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentSet;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.GeneTreeRefinement;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.IOException;
import java.io.StringReader;
import java.util.*;
import java.util.prefs.BackingStoreException;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/26/11
 * Time: 3:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProcessGT extends CommandBaseFileOut
{
    LinkedList<NetworkNonEmpty> _speciesTrees = new LinkedList<NetworkNonEmpty>();

    LinkedList<NetworkNonEmpty> _geneTrees = new LinkedList<NetworkNonEmpty>();

    private  double _bootstrap = 1;

    private HashMap<String, List<String>> _taxonMap = null;

    private boolean _geneTreesRooted = true;

    ProcessGT(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
              Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    protected int getMinNumParams()
    {
        return 2;
    }

    protected int getMaxNumParams()
    {
        return 8;
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;


        ParameterIdentSet speciesTreesIdents = super.assertParameterIdentSet(0);
        if(speciesTreesIdents != null)
        {
            _speciesTrees = assertNetworksExist(speciesTreesIdents);
            noError = _speciesTrees != null && noError;
        }

        ParameterIdentSet geneTreesIdents = super.assertParameterIdentSet(1);
        if(geneTreesIdents != null)
        {
            _geneTrees = assertNetworksExist(geneTreesIdents);
            noError = _geneTrees != null && noError;
        }

         ParamExtractor uParam = new ParamExtractor("u", this.params, this.errorDetected);
        if(uParam.ContainsSwitch)
        {
            _geneTreesRooted = false;
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


        ParamExtractorAllelListMap aParam = new ParamExtractorAllelListMap("a", this.params, this.errorDetected);

        if(aParam.ContainsSwitch)
        {
            noError = noError && aParam.IsValidMap;

            if(aParam.IsValidMap)
            {
                _taxonMap = aParam.ValueMap;
            }
        }

        checkAndSetOutFile(bParam, aParam);

        return noError;
    }

    @Override
    protected String produceResult()
    {
        List<Tree> speciesTrees = null;
		List<Tree> geneTrees = new ArrayList<Tree>();
		Map<String,List<String>> taxonMap = null;

        for(NetworkNonEmpty geneTree : _geneTrees)
        {
            String eNewickString = NetworkTransformer.toENewick(geneTree);
		    NewickReader nr = new NewickReader(new StringReader(eNewickString));
            STITree<Double> gt = new STITree<Double>(true);
            try
            {
			    nr.readTree(gt);
            }
            catch(Exception e)
            {
                throw new RuntimeException(e);
            }
            geneTrees.add(gt);
		}

        speciesTrees = new ArrayList<Tree>();
        for(NetworkNonEmpty speciesTree : _speciesTrees)
        {
            String eNewickString = NetworkTransformer.toENewick(speciesTree);
            NewickReader nr = new NewickReader(new StringReader(eNewickString));
            try
            {
			    speciesTrees.add(nr.readTree());
            }
            catch(Exception e)
            {
                throw new RuntimeException(e);
            }
		}

        StringBuffer result = new StringBuffer();
        int index = 1;
		for(Tree st: speciesTrees){
			GeneTreeRefinement.processGeneTrees(geneTrees, st, taxonMap, _geneTreesRooted, _bootstrap, null);
            String stString = st.toStringWD();
            this.treeGenerated(stString);
			result.append("\nSpecies_Tree#" + (index++ ) + " = " + stString + "\n");
			result.append("Resulting gene trees:");
			for(Tree gt: geneTrees){
                String geneTreeString = gt.toString();
				result.append("\n" + geneTreeString);
                this.treeGenerated(geneTreeString);
			}
		}

        return result.toString();


    }

}
