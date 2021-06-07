package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityIntegrated;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by hunter on 6/2/18.
 */
@CommandName("test_NCM")
public class Test_NCM extends CommandBaseFileOut {
    Double _priorBlMean = 1.;
    Double _priorBlMax = Double.POSITIVE_INFINITY;
    List<Network> _networks = new ArrayList<>();
    List<Tree> _trees = new ArrayList<>();

    public Test_NCM(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected String produceResult() {
        GeneTreeProbability gtp = new GeneTreeProbability();
        GeneTreeProbabilityIntegrated gtpi = new GeneTreeProbabilityIntegrated();
        gtpi.setBranchLengthExponentialPrior(_priorBlMean, _priorBlMax);
        String outString = "\n\n";
        for (Network n: _networks) {
            double logLikelihood = gtpi.getTotalLogLikelihood(n, _trees);
            double logLikelihoodYun = 0;
            List<Double> ls = gtp.calculateGTDistribution(n, _trees, null, false);
            for (Double likelihood: ls) {
                logLikelihoodYun += Math.log(likelihood);
            }
            outString += "network " + n.toString() + "\nNCM log likelihood = " + logLikelihood + "\nYun log likelihood = " + logLikelihoodYun + "\n\n";
        }
        return outString;

    }

    @Override
    protected int getMinNumParams() {
        return 2;
    }

    @Override
    protected int getMaxNumParams() {
        return 4;
    }

    @Override
    protected boolean checkParamsForCommand() {
        boolean noError = true;

        ParameterIdentList networkParam = this.assertParameterIdentList(0);
        for (String ident: networkParam.Elements) {
            NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
            _networks.add(transformer.makeNetwork(this.sourceIdentToNetwork.get(ident)));
        }

        ParameterIdentList treesParam = this.assertParameterIdentList(1);
        for (String ident: treesParam.Elements) {
            NetworkNonEmpty geneTree = this.sourceIdentToNetwork.get(ident);
            String phylonetGeneTree = NetworkTransformer.toENewickTree(geneTree);
            NewickReader nr = new NewickReader(new StringReader(phylonetGeneTree));
            STITree<Double> newtr = new STITree<Double>(true);
            // todo: check the tree is rooted
            try {
                nr.readTree(newtr);
            } catch (Exception e) {
                errorDetected.execute(e.getMessage(),
                        this._motivatingCommand.getLine(), this._motivatingCommand.getColumn());
            }
            Trees.removeBinaryNodes(newtr);
            for (TNode node : newtr.getNodes()) {
                node.setParentDistance(TNode.NO_DISTANCE);
            }
            _trees.add(newtr);
        }

        ParamExtractor meanParam = new ParamExtractor("mean", this.params, this.errorDetected);
        if (meanParam.ContainsSwitch) {
            _priorBlMean = Double.parseDouble(meanParam.PostSwitchValue);
        }

        ParamExtractor maxParam = new ParamExtractor("max", this.params, this.errorDetected);
        if (maxParam.ContainsSwitch) {
            _priorBlMax = Double.parseDouble(maxParam.PostSwitchValue);
        }
        noError = noError & checkForUnknownSwitches("mean", "max");
        checkAndSetOutFile(meanParam, maxParam);
        return noError;
    }
}
