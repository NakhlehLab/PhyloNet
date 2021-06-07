package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.Program;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.summary.SummaryBL;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 9/21/17
 * Time: 3:10 PM
 * To change this template use File | Settings | File Templates.
 */

@CommandName("SummarizeMCMCResults")
public class SummarizeMCMCResults extends CommandBaseFileOut {

    private String _trueNetwork = null;
    private String _outputFile = null;
    private String _pyFile = null;
    private String _videoFile = null;
    private Long _chainlen = null;
    private Long _burnin = null;
    private Long _sample_freq = null;
    private List<String> _infiles = null;
    private List<String> _leafOrder = new ArrayList<>();
    private String _mode = null;


    public SummarizeMCMCResults(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                 Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
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

        // script file
        ParamExtractor pyfileParam = new ParamExtractor("pyfile", this.params, this.errorDetected);
        if(pyfileParam.ContainsSwitch){
            if(pyfileParam.PostSwitchParam != null) {
                if (noError) {
                    _pyFile = pyfileParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -pyfile.",
                        pyfileParam.SwitchParam.getLine(), pyfileParam.SwitchParam.getColumn());
            }
        }

        // video file
        ParamExtractor videofileParam = new ParamExtractor("videofile", this.params, this.errorDetected);
        if(videofileParam.ContainsSwitch){
            if(videofileParam.PostSwitchParam != null) {
                if (noError) {
                    _videoFile = videofileParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -videofile.",
                        videofileParam.SwitchParam.getLine(), videofileParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor dataParam = new ParamExtractor("taxa", this.params, this.errorDetected);
        if(dataParam.ContainsSwitch){
            if(dataParam.PostSwitchParam != null) {
                try {
                    if(!(dataParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    ParameterIdentList taxa = (ParameterIdentList) dataParam.PostSwitchParam;
                    for(String taxon: taxa.Elements){
                        if (noError) {
                            _leafOrder.add(taxon);
                        }
                    }
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized taxa " + dataParam.PostSwitchValue,
                            dataParam.PostSwitchParam.getLine(), dataParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -taxa." + dataParam.PostSwitchValue,
                        dataParam.PostSwitchParam.getLine(), dataParam.PostSwitchParam.getColumn());
            }
        }

        // output file
        ParamExtractor outfileParam = new ParamExtractor("outfile", this.params, this.errorDetected);
        if(outfileParam.ContainsSwitch){
            if(outfileParam.PostSwitchParam != null) {
                if (noError) {
                    _outputFile = outfileParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -outfile.",
                        outfileParam.SwitchParam.getLine(), outfileParam.SwitchParam.getColumn());
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

        // chain length
        ParamExtractor clParam = new ParamExtractor("cl", this.params, this.errorDetected);
        if(clParam.ContainsSwitch){
            if(clParam.PostSwitchParam != null) {
                try {
                    _chainlen = Long.parseLong(clParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized chain length " + clParam.PostSwitchValue,
                            clParam.PostSwitchParam.getLine(), clParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -cl.",
                        clParam.SwitchParam.getLine(), clParam.SwitchParam.getColumn());
            }
        }

        // burn-in length
        ParamExtractor blParam = new ParamExtractor("bl", this.params, this.errorDetected);
        if(blParam.ContainsSwitch){
            if(blParam.PostSwitchParam != null) {
                try  {
                    _burnin = Long.parseLong(blParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized burnin length " + blParam.PostSwitchValue,
                            blParam.PostSwitchParam.getLine(), blParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -bl.",
                        blParam.SwitchParam.getLine(), blParam.SwitchParam.getColumn());
            }
        }

        // sample frequency
        ParamExtractor sfParam = new ParamExtractor("sf", this.params, this.errorDetected);
        if(sfParam.ContainsSwitch){
            if(sfParam.PostSwitchParam != null) {
                try {
                    _sample_freq = Long.parseLong(sfParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized sample frequency " + sfParam.PostSwitchValue,
                            sfParam.PostSwitchParam.getLine(), sfParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sf.",
                        sfParam.SwitchParam.getLine(), sfParam.SwitchParam.getColumn());
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
                "outfile", "infiles","pyfile","videofile",
                "truenet","taxa",
                "cl", "sf", "bl",
                "mode"
        );
        checkAndSetOutFile(
                outfileParam, tnParam, infilesParam,pyfileParam,videofileParam,
                clParam, blParam, sfParam,dataParam,
                modeParam
        );

        return  noError;
    }

    @Override
    protected String produceResult() {
        System.out.println("");

        SummaryBL sbl = new SummaryBL(_trueNetwork);
        int start = (int)(_burnin / _sample_freq) + 1;
        int end = (int)(_chainlen / _sample_freq);

        for(String file : _infiles) {
            sbl.addFile(file, true, start, end);
        }

        switch (_mode) {
            case "Topology":
                sbl.reportForTopology(_outputFile);
                break;
            case "Reticulation":
                sbl.reportForReticulation(_outputFile);
                break;
            case "Height":
                sbl.reportForHeight(_outputFile);
                break;
            case "DensityPlot":
                sbl.reportForDensityPlot(_outputFile, 1, 1);
                break;
            case "Tracer":
                sbl.reportForTracer(_outputFile, _sample_freq);
                //sbl.report( 1, 1);
                break;
            case "Video":
                sbl.reportForSlantedVideo(_pyFile, _videoFile, _leafOrder);
            default:
        }

        return "Done";
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
