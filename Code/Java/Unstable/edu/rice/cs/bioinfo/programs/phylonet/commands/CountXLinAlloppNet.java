package edu.rice.cs.bioinfo.programs.phylonet.commands;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.MDCOnAllopolyploidNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.*;

/**
 * Created by zhiyan on 7/5/20.
 */
@CommandName("deepcoalcount_alloppnet")
public class CountXLinAlloppNet extends CommandBaseFileOut{
    private HashMap<String,String> _allele2species = null;
    private HashMap<String,List<String>> _species2alleles = null;
    private NetworkNonEmpty _speciesNetwork;
    private List<NetworkNonEmpty> _geneTrees;
    private double _bootstrap = 100;

    public CountXLinAlloppNet(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                              Map<String,NetworkNonEmpty>  sourceIdentToNetwork,
                              Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams(){
        return 2;
    }

    @Override
    protected int getMaxNumParams(){
        return 9;
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

            ParamExtractor aParam = new ParamExtractor("a", this.params, this.errorDetected);
            if(aParam.ContainsSwitch){
                ParamExtractorAllelListMap alParam = new ParamExtractorAllelListMap("a", this.params, this.errorDetected);
                noError = noError && alParam.IsValidMap;
                if(alParam.IsValidMap){
                    _species2alleles = alParam.ValueMap;
                }

                ParamExtractorAllelMap aaParam = new ParamExtractorAllelMap("a", this.params, this.errorDetected);
                if(aaParam.IsValidMap){
                    _allele2species = aaParam.ValueMap;
                }
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

            noError = noError && checkForUnknownSwitches("a","b");
            checkAndSetOutFile(aParam);
        }



        return  noError;
    }

    @Override
    protected String produceResult() {
        StringBuffer result = new StringBuffer();

        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        for(NetworkNonEmpty geneTree : _geneTrees){
            double weight = geneTree.TreeProbability.execute(new TreeProbabilityAlgo<Double, RuntimeException>() {
                @Override
                public Double forEmpty(TreeProbabilityEmpty empty) {
                    return 1.0;
                }

                @Override
                public Double forNonEmpty(TreeProbabilityNonEmpty nonEmpty) {
                    return Double.parseDouble(nonEmpty.ProbString);
                }
            });

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
            Trees.removeBinaryNodes(newtr);
            if(_bootstrap<100){
                if(Trees.handleBootStrapInTree(newtr, _bootstrap)==-1){
                    throw new IllegalArgumentException("Input gene tree " + newtr + " have nodes that don't have bootstrap value");
                }

            }
            for(TNode node: newtr.getNodes()){
                node.setParentDistance(TNode.NO_DISTANCE);
            }

            String exp = Trees.getLexicographicNewickString(newtr, _allele2species);
            MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
            if(existingTuple==null){
                existingTuple = new MutableTuple<Tree, Double>(newtr, weight);
                exp2tree.put(exp, existingTuple);
            }
            else{
                existingTuple.Item2 += weight;
            }
        }

        List<MutableTuple<Tree,Double>> geneTrees = new ArrayList<MutableTuple<Tree, Double>>();
        geneTrees.addAll(exp2tree.values());

        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        Network speciesNetwork = transformer.makeNetwork(_speciesNetwork);

        int[] xlArray = new int[geneTrees.size()];

        MDCOnAllopolyploidNetwork mdc = new MDCOnAllopolyploidNetwork();
        mdc.countExtraCoal(speciesNetwork, geneTrees, _allele2species, xlArray);

        double total = 0;
        Iterator<MutableTuple<Tree,Double>> gtIt = geneTrees.iterator();
        for(int xl: xlArray){
            MutableTuple<Tree,Double> gt = gtIt.next();
            total += xl * gt.Item2;
            //result.append("\n[x" + gt.Item2 + "] " + gt.Item1.toString() + ": " + xl);
        }
        result.append("\n" + "Total number of extra lineages: " + total);

        return result.toString();
    }
}
