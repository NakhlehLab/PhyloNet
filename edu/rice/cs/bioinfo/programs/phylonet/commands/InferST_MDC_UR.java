package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCURInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCURInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/10/11
 * Time: 1:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferST_MDC_UR extends InferSTBase
{
    private Map<String,String> _taxonMap;
    private boolean _explore = false;
    private double _proportion = 0;
	private boolean _exhaust = false;
	private double _bootstrap = 1;
	private double _time = -1;
	private boolean _unresolved = false;

    InferST_MDC_UR(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 12;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand()
    {
       boolean noError = super.checkParamsForCommand();

        ParamExtractor eParam = new ParamExtractor("e", this.params, this.errorDetected);
        if(eParam.ContainsSwitch)
        {
            if(eParam.PostSwitchParam != null)
            {
                try
                {
                    _proportion = Double.parseDouble(eParam.PostSwitchValue);
                    _explore = true;
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unknown proportion value '" + eParam.PostSwitchValue + "'.", eParam.PostSwitchParam.getLine(), eParam.PostSwitchParam.getColumn());
                    noError = false;
                }
            }
            else
            {
                errorDetected.execute("Expected value after proportion switch.", eParam.SwitchParam.getLine(), eParam.SwitchParam.getColumn());
                noError = false;
            }
        }

        ParamExtractor tParam = new ParamExtractor("t", this.params, this.errorDetected);
        if(tParam.ContainsSwitch)
        {
            if(tParam.PostSwitchParam != null)
            {
                try
                {
                    _time = Double.parseDouble(tParam.PostSwitchValue);
                }
                catch(NumberFormatException e)
                {
                    errorDetected.execute("Unknown time value '" + tParam.PostSwitchValue + "'.", tParam.PostSwitchParam.getLine(), tParam.PostSwitchParam.getColumn());
                    noError = false;
                }
            }
            else
            {
                errorDetected.execute("Expected value after time switch.", tParam.SwitchParam.getLine(), tParam.SwitchParam.getColumn());
                noError = false;
            }
        }

        ParamExtractor xParam = new ParamExtractor("x", this.params, this.errorDetected);
        if(xParam.ContainsSwitch)
        {
            _exhaust = true;
        }

        ParamExtractor urParam = new ParamExtractor("ur", this.params, this.errorDetected);
        if(urParam.ContainsSwitch)
        {
            _unresolved = true;
        }



        ThresholdResult tr = this.assignThreshold(1);
        _bootstrap = tr.Threshold;
        noError = noError && tr.NoError;

        TaxonMapResult result = assignTaxonMap();
        noError = noError && result.NoError;
        _taxonMap = result.TaxonMap;

        this.checkAndSetOutFile(eParam, tr.Extractor, result.Extractor, tParam);


       return  noError;
    }

    @Override
    protected String produceResult()
    {
        if(_geneTrees == null)
        {
            throw new IllegalStateException();
        }

        StringBuffer result = new StringBuffer();

        List<Tree> trees = GetGeneTreesAsTreeList();

        MDCURInference_DP inference = new MDCURInference_DP();

        List<Solution> solutions;
		if(_taxonMap == null){
			solutions = inference.inferSpeciesTree(trees,_explore,_proportion,_exhaust,_bootstrap,_unresolved,_time);
		}
		else
            solutions = inference.inferSpeciesTree(trees,_taxonMap,_explore,_proportion,_exhaust,_bootstrap,_unresolved,_time);

			for(Solution s: solutions){
				result.append("\n" + s._st.toString()+" "+s._totalCoals+" extra lineages in total");
			}


        return result.toString();
    }
}
