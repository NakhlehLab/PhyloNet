package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 5/9/17
 * Time: 11:15 AM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("SimBiMarkersinNetwork")
public class SimBiMarkersinNetwork extends CommandBaseFileOut{

    private boolean _diploid = false;
    private Integer _polyploid = null; // TODO: not fully tested
    private Character _dominant = null;
    private boolean _useOnlyPolymorphic = false;
    private Long _seed = null;
    private String _trueNetwork = null;
    private Integer _n = null;
    private Integer _nSitesPerGT = null;
    private Boolean _rateVarAcrossMarkers = null;
    private Boolean _rateVarAcrossLineages = null;
    private Double _invariantSitesProb = null;
    private Double _gammaShape = null;
    private Integer _gammaCat = null;
    private Map<String, List<String>> _species2alleles = null;
    private double _pi0 = 0.5;
    private String _outFilePath = null;
    private Double _theta = null;


    public SimBiMarkersinNetwork(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                          Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 2;
    }

    @Override
    protected int getMaxNumParams() {
        return 50;
    }

    @Override
    protected boolean checkParamsForCommand() {
        boolean noError = true;

        //diploid
        ParamExtractor diploidParam = new ParamExtractor("diploid", this.params, this.errorDetected);
        if(diploidParam.ContainsSwitch) {
            _diploid = true;
        }

        //dominant
        ParamExtractor dominantParam = new ParamExtractor("dominant", this.params, this.errorDetected);
        if(dominantParam.ContainsSwitch){
            _dominant = '0';
            if(!_diploid) {
                errorDetected.execute("Only diploid has dominant markers.",
                        dominantParam.PostSwitchParam.getLine(), dominantParam.PostSwitchParam.getColumn());
            }
        }

        // polyploid
        ParamExtractor polyParam = new ParamExtractor("polyploid", this.params, this.errorDetected);
        if(polyParam.ContainsSwitch){
            if(polyParam.PostSwitchParam != null) {
                try {
                    _polyploid = Integer.parseInt(polyParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized polyploid " + polyParam.PostSwitchValue,
                            polyParam.PostSwitchParam.getLine(), polyParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -polyploid.",
                        polyParam.SwitchParam.getLine(), polyParam.SwitchParam.getColumn());
            }
        }

        // only polymorphic
        ParamExtractor onlyPolyParam = new ParamExtractor("op", this.params, this.errorDetected);
        if(onlyPolyParam.ContainsSwitch){
            _useOnlyPolymorphic = true;
        }

        // seed
        ParamExtractor sdParam = new ParamExtractor("sd", this.params, this.errorDetected);
        if(sdParam.ContainsSwitch){
            if(sdParam.PostSwitchParam != null) {
                try {
                    _seed = Long.parseLong(sdParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized random seed " + sdParam.PostSwitchValue,
                            sdParam.PostSwitchParam.getLine(), sdParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sd.",
                        sdParam.SwitchParam.getLine(), sdParam.SwitchParam.getColumn());
            }
        }

        // using coalescent unit and give theta
        ParamExtractor cuParam = new ParamExtractor("cu", this.params, this.errorDetected);
        if(cuParam.ContainsSwitch){
            if(cuParam.PostSwitchParam != null) {
                try {
                    _theta = Double.parseDouble(cuParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized theta " + cuParam.PostSwitchValue,
                            cuParam.PostSwitchParam.getLine(), cuParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -cu.",
                        cuParam.SwitchParam.getLine(), cuParam.SwitchParam.getColumn());
            }
        }

        // invariant sites
        ParamExtractor invariantParam = new ParamExtractor("i", this.params, this.errorDetected);
        if(invariantParam.ContainsSwitch){
            if(invariantParam.PostSwitchParam != null) {
                try {
                    _invariantSitesProb = Double.parseDouble(invariantParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized invariant sites prob " + invariantParam.PostSwitchValue,
                            invariantParam.PostSwitchParam.getLine(), invariantParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -i.",
                        invariantParam.SwitchParam.getLine(), invariantParam.SwitchParam.getColumn());
            }
        }

        // gamma categories
        ParamExtractor gammaCatParam = new ParamExtractor("g", this.params, this.errorDetected);
        if(gammaCatParam.ContainsSwitch){
            if(gammaCatParam.PostSwitchParam != null) {
                try {
                    _gammaCat = Integer.parseInt(gammaCatParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized number of categories " + gammaCatParam.PostSwitchValue,
                            gammaCatParam.PostSwitchParam.getLine(), gammaCatParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -g.",
                        gammaCatParam.SwitchParam.getLine(), gammaCatParam.SwitchParam.getColumn());
            }
        }

        // gamma shape
        ParamExtractor gammaShapeParam = new ParamExtractor("a", this.params, this.errorDetected);
        if(gammaShapeParam.ContainsSwitch){
            if(gammaShapeParam.PostSwitchParam != null) {
                try {
                    _gammaShape = Double.parseDouble(gammaShapeParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized alpha " + gammaShapeParam.PostSwitchValue,
                            gammaShapeParam.PostSwitchParam.getLine(), gammaShapeParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -a.",
                        gammaShapeParam.SwitchParam.getLine(), gammaShapeParam.SwitchParam.getColumn());
            }
        }

        // rvl
        ParamExtractor rvlParam = new ParamExtractor("rvl", this.params, this.errorDetected);
        if(rvlParam.ContainsSwitch){
            _rateVarAcrossLineages = true;
        }

        // rvm
        ParamExtractor rvmParam = new ParamExtractor("rvm", this.params, this.errorDetected);
        if(rvmParam.ContainsSwitch){
            _rateVarAcrossMarkers = true;
        }

        // sitespergt
        ParamExtractor sitespergtParam = new ParamExtractor("sitespergt", this.params, this.errorDetected);
        if(sitespergtParam.ContainsSwitch){
            if(sitespergtParam.PostSwitchParam != null) {
                try {
                    _nSitesPerGT = Integer.parseInt(sitespergtParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized number of sites " + sitespergtParam.PostSwitchValue,
                            sitespergtParam.PostSwitchParam.getLine(), sitespergtParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sitespergt.",
                        sitespergtParam.SwitchParam.getLine(), sitespergtParam.SwitchParam.getColumn());
            }
        }

        // n
        ParamExtractor numParam = new ParamExtractor("num", this.params, this.errorDetected);
        if(numParam.ContainsSwitch){
            if(numParam.PostSwitchParam != null) {
                try {
                    _n = Integer.parseInt(numParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized number of sites " + numParam.PostSwitchValue,
                            numParam.PostSwitchParam.getLine(), numParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -num.",
                        numParam.SwitchParam.getLine(), numParam.SwitchParam.getColumn());
            }
        }

        // pi0
        ParamExtractor pi0Param = new ParamExtractor("pi0", this.params, this.errorDetected);
        if(pi0Param.ContainsSwitch){
            if(pi0Param.PostSwitchParam != null) {
                try {
                    _pi0 = Double.parseDouble(pi0Param.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized stationary probabiltiy " + pi0Param.PostSwitchValue,
                            pi0Param.PostSwitchParam.getLine(), pi0Param.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -pi0.",
                        pi0Param.SwitchParam.getLine(), pi0Param.SwitchParam.getColumn());
            }
        }


        // taxon map
        ParamExtractor tmParam = new ParamExtractor("tm", this.params, this.errorDetected);
        if(tmParam.ContainsSwitch){
            ParamExtractorAllelListMap aaParam = new ParamExtractorAllelListMap("tm", this.params, this.errorDetected);
            noError = noError && aaParam.IsValidMap;
            if(aaParam.IsValidMap){
                _species2alleles = aaParam.ValueMap;
            }
        }

        // true network
        ParamExtractor snParam = new ParamExtractor("truenet", this.params, this.errorDetected);
        if(snParam.ContainsSwitch){
            if(snParam.PostSwitchParam != null) {
                if (noError) {
                    _trueNetwork = snParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -truenet.",
                        snParam.SwitchParam.getLine(), snParam.SwitchParam.getColumn());
            }
        }

        // output
        ParamExtractor fileParam = new ParamExtractor("out", this.params, this.errorDetected);
        if(fileParam.ContainsSwitch){
            if(fileParam.PostSwitchParam != null) {
                if (noError) {
                    _outFilePath = fileParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -out.",
                        fileParam.SwitchParam.getLine(), fileParam.SwitchParam.getColumn());
            }
        }

        noError = noError && checkForUnknownSwitches(
                "diploid",
                "dominant","polyploid","pi0","cu",
                "op", "sd", "num", "tm", "truenet", "out",
                "sitespergt", "rvl", "rvm",
                "a", "i", "g"
        );
        checkAndSetOutFile(
                diploidParam,dominantParam,polyParam,onlyPolyParam,pi0Param,cuParam,
                numParam, sdParam, tmParam,snParam, fileParam,
                sitespergtParam, rvlParam, rvmParam,
                gammaShapeParam, gammaCatParam, invariantParam
        );

        return  noError;

    }

    @Override
    protected String produceResult() {
        double pi0 = _pi0;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = _useOnlyPolymorphic;
        int numSites = _n;

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, _seed);
        simulator._diploid = _diploid;
        simulator._dominant = _dominant != null;
        simulator._polyploid = _polyploid;

        if(_nSitesPerGT != null) {
            simulator.setSitesPerGT(_nSitesPerGT);
        }

        if(_rateVarAcrossLineages != null) {
            simulator.setRateVarAcrossLineages(_rateVarAcrossLineages);
        }

        if(_rateVarAcrossMarkers != null) {
            simulator.setRateVarAcrossMarkers(_rateVarAcrossMarkers);
        }

        if(_invariantSitesProb != null) {
            simulator.setInvariantSitesProb(_invariantSitesProb);
        }

        if(_gammaCat != null) {
            simulator.setDiscreteGamma(_gammaShape, _gammaCat, _seed);
        }

        String s = _trueNetwork;
        double rootPopSize = 0.0;
        if(s.startsWith("[")) {
            rootPopSize = Double.parseDouble(s.substring(1, s.indexOf("]")));
            s = s.substring(s.indexOf("]") + 1);
        }
        Network<Object> trueNetwork = edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(s);
        trueNetwork.getRoot().setRootPopSize(rootPopSize);

        if(_theta != null) {
            for (Object nodeObject : edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.postTraversal(trueNetwork)) {
                NetNode node = (NetNode) nodeObject;
                for (Object parentObject : node.getParents()) {
                    NetNode parent = (NetNode) parentObject;
                    node.setParentSupport(parent, _theta);
                    node.setParentDistance(parent, node.getParentDistance(parent) * _theta / 2.0);
                }
            }
            trueNetwork.getRoot().setRootPopSize(_theta);
        }

        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, _species2alleles, numSites, !useOnlyPolymorphic);

        StringBuilder result = new StringBuilder();
        result.append("\n\n");

        result.append("True Network: " + _trueNetwork + "\n");
        result.append("Seed: " + _seed + "\n");
        result.append("Stationary Probability: " + pi0 + " " + pi1 + "\n");
        result.append("Only polymorphic site: " + useOnlyPolymorphic + "\n");
        result.append("Diploid: " + _diploid + "\n");
        result.append("Dominant marker: " + _dominant + "\n");

        result.append("\n");

        List<String> taxa = new ArrayList<>(onesnp.keySet());
        Collections.sort(taxa);

        if(_outFilePath == null) {

            for (String allele : taxa) {
                result.append(allele);
                result.append(" ");
                result.append(onesnp.get(allele));
                result.append("\n");
            }
        } else {
            result.append("Write to: ");
            result.append(_outFilePath);
            result.append("\n");
            try {
                PrintWriter out = new PrintWriter(_outFilePath);
                for (String taxon : taxa) {
                    out.println(taxon + " " + onesnp.get(taxon));
                }
                out.close();
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
        }

        return result.toString();
    }
}
