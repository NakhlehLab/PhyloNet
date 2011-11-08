package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentSet;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.DeepCoalescencesCounter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.StringReader;
import java.security.cert.TrustAnchor;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 4:45 PM
 * To change this template use File | Settings | File Templates.
 */
public class DeepCoalCount extends CommandBaseFileOut {

    private Object taxonMap = null;

    private boolean  _treatGeneTreesAsRooted = true;

    private double _bootstrap = 1.0;

    private List<NetworkNonEmpty> _speciesTrees;

    private List<NetworkNonEmpty> _geneTrees;

    DeepCoalCount(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    protected int getMinNumParams()
    {
        return 2;
    }

    protected int getMaxNumParams()
    {
        return 6;
    }

    public boolean checkParamsForCommand() {

        boolean noError = true;

        ParameterIdentSet speciesTreeParam = this.assertParameterIdentSet(0);
        ParameterIdentSet geneTreeParam    = this.assertParameterIdentSet(1);

        for(int i = 2; i<this.params.size(); i++)
        {
            Parameter ithParam = params.get(i);
            String value = params.get(i).execute(GetSimpleParamValue.Singleton, null);

            if(value != null)
            {
                if(value.toLowerCase().equals("-b"))
                {
                    if(i + 1 < params.size())
                    {
                        Parameter ithPlusOneParam = params.get(i + 1);

                    }
                    else
                    {
                        errorDetected.execute("Expected threshold value after -b switch.", ithParam.getLine(), ithParam.getColumn());
                    }
                }
            }
            else
            {
                errorDetected.execute("Expected -b or -a switches.", ithParam.getLine(), ithParam.getColumn());
            }

        }


        if(noError)
        {
            return checkContext(speciesTreeParam, geneTreeParam);
        }
        else
        {
            return noError;
        }

    }

    private boolean checkContext(ParameterIdentSet speciesTreeParam, ParameterIdentSet geneTreeParam ) {

        boolean noError = true;

        _speciesTrees = new LinkedList<NetworkNonEmpty>();
        for(String ident : speciesTreeParam.Elements)
        {
            noError = noError && this.assertNetworkExists(ident, speciesTreeParam.getLine(), speciesTreeParam.getColumn());
            if(noError)
            {
                _speciesTrees.add(this.sourceIdentToNetwork.get(ident));
            }
        }

        _geneTrees = new LinkedList<NetworkNonEmpty>();
        for(String ident : geneTreeParam.Elements)
        {
            noError = noError && this.assertNetworkExists(ident, geneTreeParam.getLine(), geneTreeParam.getColumn());
            if(noError)
            {
                _geneTrees.add(this.sourceIdentToNetwork.get(ident));
            }
        }

        if(!noError)
        {
            _speciesTrees = null;
            _geneTrees = null;

        }

        return noError;

    }

    @Override
    protected String produceResult() {

        StringBuffer result = new StringBuffer();

       if(taxonMap == null)
       {
           List<Tree> geneTrees = new LinkedList<Tree>();

           for(NetworkNonEmpty geneTree : _geneTrees)
           {
               String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
               NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
			   STITree<Double> gt = new STITree<Double>(true);
			   try
               {
                   nr.readTree(gt);
                   geneTrees.add(gt);
               }
               catch(Exception e)
               {
                    errorDetected.execute(e.getMessage(), -1, -1);
               }

           }

           int index = 1;
           for(NetworkNonEmpty st : _speciesTrees)
           {
               String phylonetSpeciesTree = NetworkTransformer.toENewickTree(st);
               NewickReader nr = new NewickReader(new StringReader(phylonetSpeciesTree));

               try
               {
                    Tree speciesTree = nr.readTree();
                    int coalNum = DeepCoalescencesCounter.countExtraCoal(geneTrees, speciesTree, _treatGeneTreesAsRooted, _bootstrap);
                    result.append("\nSpecies_Tree#" + (index++ ) + " = " + speciesTree.toStringWD() + "\n");
			        result.append("Total number of extra lineages: " + coalNum);
               }
               catch(Exception e)
               {
                   errorDetected.execute(e.getMessage(), -1, -1);
               }
           }
       }

        return result.toString();

    }
}
