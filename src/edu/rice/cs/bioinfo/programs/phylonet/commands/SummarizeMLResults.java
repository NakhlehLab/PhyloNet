package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.Program;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 1/3/18
 * Time: 9:06 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("SummarizeMLResults")
public class SummarizeMLResults extends CommandBaseFileOut {

    private List<String> _infiles = null;
    private String _mode = null;
    private String _trueNetwork = null;



    public SummarizeMLResults(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                Proc3<String, Integer, Integer> errorDetected, RichNewickReader<edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 2;
    }

    @Override
    protected int getMaxNumParams() {
        return 30;
    }

    @Override
    protected boolean checkParamsForCommand() {
        boolean noError = true;

        _infiles = parseInputFiles(Program.inputNexusFileName);

        // mode
        ParamExtractor modeParam = new ParamExtractor("mode", this.params, this.errorDetected);
        if(modeParam.ContainsSwitch){
            if(modeParam.PostSwitchParam != null) {
                if (noError) {
                    _mode = modeParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -mode.",
                        modeParam.SwitchParam.getLine(), modeParam.SwitchParam.getColumn());
            }
        }

        // true network
        ParamExtractor tnParam = new ParamExtractor("truenet", this.params, this.errorDetected);
        if(tnParam.ContainsSwitch){
            if(tnParam.PostSwitchParam != null) {
                if (noError) {
                    _trueNetwork = tnParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -truenet.",
                        tnParam.SwitchParam.getLine(), tnParam.SwitchParam.getColumn());
            }
        }

        // input files
        ParamExtractor infilesParam = new ParamExtractor("infiles", this.params, this.errorDetected);
        if(infilesParam.ContainsSwitch) {
            if(infilesParam.PostSwitchParam != null) {
                try {
                    if(!(infilesParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    _infiles = new ArrayList<>();
                    ParameterIdentList temps = (ParameterIdentList) infilesParam.PostSwitchParam;
                    for(String tExp: temps.Elements){
                        _infiles.add(tExp.trim());
                    }
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -infiles.",
                            infilesParam.PostSwitchParam.getLine(), infilesParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -infiles.",
                        infilesParam.SwitchParam.getLine(), infilesParam.SwitchParam.getColumn());
            }
        }

        noError = noError && checkForUnknownSwitches(
                "truenet", "infiles",
                "mode"
        );
        checkAndSetOutFile(
                tnParam, infilesParam,
                modeParam
        );

        return noError;
    }

    @Override
    protected String produceResult() {
        System.out.println("");

        switch (_mode) {
            case "Report":
                report();
                break;
            default:
        }

        return "Done";
    }

    private void report() {
        Network trueNetwork = Networks.readNetworkWithRootPop(_trueNetwork);
        int topOptimalCorrect = 0;
        int topFiveOptimalCorrect = 0;

        for(String infile : _infiles) {
            List<String> optimalNetworkStrings = new ArrayList<>();
            int finalRound = 0;

            try {
                BufferedReader in = new BufferedReader(new FileReader(infile));

                String s;
                int numRound = -1;

                while((s = in.readLine()) != null) {
                    if(s.contains("Results after run #")) {
                        numRound = Integer.parseInt(s.substring(s.indexOf('#') + 1));
                        optimalNetworkStrings.clear();
                        continue;
                    }

                    if(numRound != -1) {
                        Pattern pattern = Pattern.compile("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?");
                        Matcher matcher = pattern.matcher(s);
                        if(matcher.find() && matcher.start() == 0) {
                            // parse optimal network lists
                            Pattern pattern2 = Pattern.compile(":");
                            Matcher matcher2 = pattern2.matcher(s);
                            if(matcher2.find() && matcher2.find()) {
                                optimalNetworkStrings.add(s.substring(matcher2.end() + 1));
                            }
                        } else if(s.startsWith("Running Time")) {
                            // end of one round
                            finalRound = numRound;
                            numRound = -1;
                        }
                    }
                }

                in.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            System.out.println(infile + " : final round: " + finalRound);

            Network topOptimalNetwork = Networks.readNetworkWithRootPop(optimalNetworkStrings.get(0));
            if(Networks.hasTheSameTopology(trueNetwork, topOptimalNetwork)) {
                System.out.println("Best network correct");
                topOptimalCorrect++;
            } else {
                for(int i = 1 ; i < 5 ; i++) {
                    Network network = Networks.readNetworkWithRootPop(optimalNetworkStrings.get(i));
                    if(Networks.hasTheSameTopology(trueNetwork, network)) {
                        System.out.println("Top 5 network correct");
                        topFiveOptimalCorrect++;
                    }
                }
            }

        }

        System.out.println("Best network correct: " + topOptimalCorrect);
        System.out.println("Top 5 network correct: " + topFiveOptimalCorrect);
    }

    private List<String> parseInputFiles(String file) {
        List<String> result = new ArrayList<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String s;
            boolean begin = false;

            while((s = br.readLine().trim()) != null) {
                if(s.toLowerCase().equals("begin sets;")) {
                    begin = true;
                    continue;
                }
                if(begin && s.toLowerCase().equals("end;")) {
                    break;
                }
                if(begin) {
                    result.add(s);
                }
            }

            br.close();
        } catch (Exception ex) {
            result = null;
        }
        return result;
    }
}
