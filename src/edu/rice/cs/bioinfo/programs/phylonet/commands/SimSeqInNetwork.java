package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSeqInNetworkByMSSeqGen;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 3/24/18
 * Time: 5:40 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("SimSeqInNetwork")
public class SimSeqInNetwork extends CommandBaseFileOut{
    private String _trueNetwork = null;
    private Integer _numGT = null;
    private Integer _nSitesPerGT = null;
    private String _outFilePath = null;
    private String _gtFilePath = null;
    private String _msPath = null;
    private String _seqgenPath = null;
    private Map<String, List<String>> _species2alleles = null;
    private List<Double> _rates = null;
    private List<Double> _bases = null;
    private Double _theta = null;


    public SimSeqInNetwork(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
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
        ParamExtractor numParam = new ParamExtractor("numgt", this.params, this.errorDetected);
        if(numParam.ContainsSwitch){
            if(numParam.PostSwitchParam != null) {
                try {
                    _numGT = Integer.parseInt(numParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized number of sites " + numParam.PostSwitchValue,
                            numParam.PostSwitchParam.getLine(), numParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -numgt.",
                        numParam.SwitchParam.getLine(), numParam.SwitchParam.getColumn());
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

        // output gt
        ParamExtractor gtfileParam = new ParamExtractor("gtfile", this.params, this.errorDetected);
        if(gtfileParam.ContainsSwitch){
            if(gtfileParam.PostSwitchParam != null) {
                if (noError) {
                    _gtFilePath = gtfileParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -gtfile.",
                        gtfileParam.SwitchParam.getLine(), gtfileParam.SwitchParam.getColumn());
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

        // taxon map
        ParamExtractor tmParam = new ParamExtractor("tm", this.params, this.errorDetected);
        if(tmParam.ContainsSwitch){
            ParamExtractorAllelListMap aaParam = new ParamExtractorAllelListMap("tm", this.params, this.errorDetected);
            noError = noError && aaParam.IsValidMap;
            if(aaParam.IsValidMap){
                _species2alleles = aaParam.ValueMap;
            }
        }

        // mspath
        ParamExtractor mspathParam = new ParamExtractor("mspath", this.params, this.errorDetected);
        if(mspathParam.ContainsSwitch){
            if(mspathParam.PostSwitchParam != null) {
                if (noError) {
                    _msPath = mspathParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -mspath.",
                        mspathParam.SwitchParam.getLine(), mspathParam.SwitchParam.getColumn());
            }
        }

        // seqgenpath
        ParamExtractor seqgenpathParam = new ParamExtractor("seqgenpath", this.params, this.errorDetected);
        if(seqgenpathParam.ContainsSwitch){
            if(seqgenpathParam.PostSwitchParam != null) {
                if (noError) {
                    _seqgenPath = seqgenpathParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -seqgenpath.",
                        seqgenpathParam.SwitchParam.getLine(), seqgenpathParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor rateParam = new ParamExtractor("rate", this.params, this.errorDetected);
        if(rateParam.ContainsSwitch) {
            if(rateParam.PostSwitchParam != null) {
                try {
                    if(!(rateParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    _rates = new ArrayList<>();
                    ParameterIdentList temps = (ParameterIdentList) rateParam.PostSwitchParam;
                    for(String tExp: temps.Elements){
                        double t = Double.parseDouble(tExp.trim());
                        if(t < 0) {
                            throw new RuntimeException();
                        }
                        _rates.add(t);
                    }
                    if(_rates.size() != 6) {
                        throw new NumberFormatException("There should be 6 rates");
                    }
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -rate.",
                            rateParam.PostSwitchParam.getLine(), rateParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -rate.",
                        rateParam.SwitchParam.getLine(), rateParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor baseParam = new ParamExtractor("base", this.params, this.errorDetected);
        if(baseParam.ContainsSwitch) {
            if(baseParam.PostSwitchParam != null) {
                try {
                    if(!(baseParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    _bases = new ArrayList<>();
                    ParameterIdentList temps = (ParameterIdentList) baseParam.PostSwitchParam;
                    for(String tExp: temps.Elements){
                        double t = Double.parseDouble(tExp.trim());
                        if(t < 0) {
                            throw new RuntimeException();
                        }
                        _bases.add(t);
                    }
                    if(_bases.size() != 4) {
                        throw new NumberFormatException("There should be 4 bases");
                    }
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -base.",
                            baseParam.PostSwitchParam.getLine(), baseParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -base.",
                        baseParam.SwitchParam.getLine(), baseParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor thetaParam = new ParamExtractor("theta", this.params, this.errorDetected);
        if(thetaParam.ContainsSwitch){
            if(thetaParam.PostSwitchParam != null) {
                try {
                    _theta = Double.parseDouble(thetaParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized theta " + thetaParam.PostSwitchValue,
                            thetaParam.PostSwitchParam.getLine(), thetaParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -theta.",
                        thetaParam.SwitchParam.getLine(), thetaParam.SwitchParam.getColumn());
            }
        }

        noError = noError && checkForUnknownSwitches(
                "numgt", "tm", "truenet", "out", "gtfile",
                "sitespergt", "mspath", "seqgenpath", "base", "rate",
                "theta"
        );
        checkAndSetOutFile(
                numParam, tmParam,snParam, fileParam, gtfileParam,
                sitespergtParam, mspathParam, seqgenpathParam, baseParam, rateParam,
                thetaParam
        );

        return  noError;
    }

    @Override
    protected String produceResult() {
        StringBuilder result = new StringBuilder();
        double base[] = new double [4];
        base[0] = _bases.get(0);
        base[1] = _bases.get(1);
        base[2] = _bases.get(2);
        base[3] = _bases.get(3);

        double rate[] = new double[6];
        rate[0] = _rates.get(0);
        rate[1] = _rates.get(1);
        rate[2] = _rates.get(2);
        rate[3] = _rates.get(3);
        rate[4] = _rates.get(4);
        rate[5] = _rates.get(5);

        Network net = edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(_trueNetwork);
        List<Map<String, String>> lociSeq = SimSeqInNetworkByMSSeqGen.generateSeqences(net, _species2alleles, _numGT, _nSitesPerGT, _theta, base, rate, _msPath, _seqgenPath, _gtFilePath);

        int i = 0;

        if(_outFilePath == null) {
            for (Map<String, String> seq : lociSeq) {
                i++;
                result.append(String.format("[loci%d, %d, ...]\n", i, seq.values().iterator().next().length()));
                for (String key : seq.keySet()) {
                    result.append(String.format("%s %s\n", key, seq.get(key)));
                }
            }
        } else {
            result.append("Write to: ");
            result.append(_outFilePath);
            result.append("\n");

            try {
                PrintWriter out = new PrintWriter(_outFilePath);
                for (Map<String, String> seq : lociSeq) {
                    i++;
                    out.print(String.format("[loci%d, %d, ...]\n", i, seq.values().iterator().next().length()));
                    for (String key : seq.keySet()) {
                        out.print(String.format("%s %s\n", key, seq.get(key)));
                    }
                }
                out.close();
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
        }
        return result.toString();
    }
}
