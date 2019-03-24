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
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SNOptions;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SNSummary;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 5/27/18
 * Time: 11:43 AM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("NetMerger")
public class NetMerger extends CommandBaseFileOut {
    private String _trueNetwork = null;
    private String _outputFile = null;
    private String _pyFile = null;
    private String _inputFolder = null;
    private Integer _chainlen = null;
    private Integer _burnin = null;
    private Integer _sample_freq = null;
    private Double _eps = 0.01;
    private String _outgroup = null;


    public NetMerger(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
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

        // mode
        ParamExtractor outgroupParam = new ParamExtractor("outgroup", this.params, this.errorDetected);
        if(outgroupParam.ContainsSwitch){
            if(outgroupParam.PostSwitchParam != null) {
                if (noError) {
                    _outgroup = outgroupParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -outgroup.",
                        outgroupParam.SwitchParam.getLine(), outgroupParam.SwitchParam.getColumn());
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
        ParamExtractor inputFolderParam = new ParamExtractor("inputfolder", this.params, this.errorDetected);
        if(inputFolderParam.ContainsSwitch){
            if(inputFolderParam.PostSwitchParam != null) {
                if (noError) {
                    _inputFolder = inputFolderParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -inputfolder.",
                        inputFolderParam.SwitchParam.getLine(), inputFolderParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor epsParam = new ParamExtractor("eps", this.params, this.errorDetected);
        if(epsParam.ContainsSwitch){
            if(epsParam.PostSwitchParam != null) {
                try {
                    _eps = Double.parseDouble(epsParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized Poisson parameter " + epsParam.PostSwitchValue,
                            epsParam.PostSwitchParam.getLine(), epsParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -eps.",
                        epsParam.SwitchParam.getLine(), epsParam.SwitchParam.getColumn());
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
                    _chainlen = Integer.parseInt(clParam.PostSwitchValue);
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
                    _burnin = Integer.parseInt(blParam.PostSwitchValue);
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
                    _sample_freq = Integer.parseInt(sfParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized sample frequency " + sfParam.PostSwitchValue,
                            sfParam.PostSwitchParam.getLine(), sfParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sf.",
                        sfParam.SwitchParam.getLine(), sfParam.SwitchParam.getColumn());
            }
        }

        noError = noError && checkForUnknownSwitches(
                "eps","inputfolder",
                "truenet",
                "cl", "sf", "bl",
                "outgroup"
        );
        checkAndSetOutFile(
                epsParam, tnParam,inputFolderParam,
                clParam, blParam, sfParam,
                outgroupParam
        );

        return  noError;
    }

    @Override
    protected String produceResult() {
        System.out.println("");

        String resultFolder = _inputFolder;
        File path = new File(resultFolder);

        List<String> filenames = new ArrayList<>();

        File [] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].toString().endsWith(".out")) {
                    //System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);
        SNOptions options = new SNOptions();
        options.outgroup = _outgroup;
        options.eps = _eps;

        SuperNetwork3.printDetails_ = true;
        SNSummary summary = Pipeline.stage2_1(filenames, _chainlen, _burnin, _sample_freq, options);

        Network inferred = summary.inferredNetwork;
        System.out.println(inferred);

        return "";
    }
}
