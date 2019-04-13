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
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.reducing.GenerateTrinets3Sets;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import org.jfree.chart.plot.PieLabelDistributor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.*;

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
    private String _nexTemplateFile = null;
    private String _gtFile = null;
    private String _tripletFile = null;
    private String _inputFolder = null;
    private Integer _chainlen = null;
    private Integer _burnin = null;
    private Integer _sample_freq = null;
    private Double _eps = 0.01;
    private String _outgroup = null;
    private String _mode = null;
    private String _order = null;
    private String _backbone = null;
    private String _cmd;
    private List<String> _species = null;
    private Map<String, List<String>> _species2alleles;


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

    String parseCmd(String file) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));

            String s;
            StringBuilder cmd = new StringBuilder();
            boolean begin = false;
            while((s = br.readLine().trim()) != null ) {
                String ss = s.replaceAll("\\s+", "").toLowerCase();
                if (ss.startsWith("beginphylonet;")) {
                    begin = true;
                }
                else if(begin) {
                    if(ss.endsWith("end;")) {
                        break;
                    }
                    cmd.append(s + " ");
                }
            }
            return cmd.toString();

        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    protected boolean checkParamsForCommand() {
        boolean noError = true;


        // taxon map
        ParamExtractor tmParam = new ParamExtractor("tm", this.params, this.errorDetected);
        if(tmParam.ContainsSwitch){
            ParamExtractorAllelListMap aaParam = new ParamExtractorAllelListMap("tm", this.params, this.errorDetected);
            noError = noError && aaParam.IsValidMap;
            if(aaParam.IsValidMap){
                _species2alleles = aaParam.ValueMap;
            }
        }



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

        // order
        ParamExtractor orderParam = new ParamExtractor("order", this.params, this.errorDetected);
        if(orderParam.ContainsSwitch){
            if(orderParam.PostSwitchParam != null) {
                if (noError) {
                    _order = orderParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -order.",
                        orderParam.SwitchParam.getLine(), orderParam.SwitchParam.getColumn());
            }
        }

        // backbone
        ParamExtractor backboneParam = new ParamExtractor("backbone", this.params, this.errorDetected);
        if(backboneParam.ContainsSwitch){
            if(backboneParam.PostSwitchParam != null) {
                if (noError) {
                    _backbone = backboneParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -backbone.",
                        backboneParam.SwitchParam.getLine(), backboneParam.SwitchParam.getColumn());
            }
        }

        // outgroup
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

        // triplet file
        ParamExtractor tripletsParam = new ParamExtractor("triplets", this.params, this.errorDetected);
        if(tripletsParam.ContainsSwitch){
            if(tripletsParam.PostSwitchParam != null) {
                if (noError) {
                    _tripletFile = tripletsParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -triplets.",
                        tripletsParam.SwitchParam.getLine(), tripletsParam.SwitchParam.getColumn());
            }
        }

        // gt file
        ParamExtractor gtsParam = new ParamExtractor("gts", this.params, this.errorDetected);
        if(gtsParam.ContainsSwitch){
            if(gtsParam.PostSwitchParam != null) {
                if (noError) {
                    _gtFile = gtsParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -gts.",
                        gtsParam.SwitchParam.getLine(), gtsParam.SwitchParam.getColumn());
            }
        }

        // nex template file
        ParamExtractor nexParam = new ParamExtractor("nex", this.params, this.errorDetected);
        if(nexParam.ContainsSwitch){
            if(nexParam.PostSwitchParam != null) {
                if (noError) {
                    _nexTemplateFile = nexParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected string after switch -nex.",
                        nexParam.SwitchParam.getLine(), nexParam.SwitchParam.getColumn());
            }
        }

        // input folder
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
                "truenet", "mode",
                "cl", "sf", "bl",
                "outgroup", "nex", "tm", "gts", "triplets",
                "backbone", "order"
        );
        checkAndSetOutFile(
                epsParam, tnParam,inputFolderParam,
                modeParam,
                clParam, blParam, sfParam,
                outgroupParam, nexParam, tmParam, gtsParam, tripletsParam,
                backboneParam, orderParam
        );

        return  noError;
    }

    private void writeTriplets() {
        if(_outgroup == null) {
            throw new RuntimeException("Outgroup is needed to generate triplets!");
        }

        // If gene trees are given, use reduced set of trinets.
        if(_gtFile != null) {
            GenerateTrinets3Sets gen = new GenerateTrinets3Sets();

            try {
                BufferedReader br = new BufferedReader(new FileReader(_gtFile));
                String s;
                while((s = br.readLine()) != null ) {
                    gen.addTree(s);
                }
                br.close();

            } catch (Exception e) {
                e.printStackTrace();
            }

            Set<String[]> sets = gen.generate();

            List<Set<String>> speciesTriplets = Pipeline.convertAlleleTripletsToSpeciesTriplets(_species2alleles, _outgroup, sets);

            try {
                PrintWriter out = new PrintWriter(_tripletFile);
                for(Set<String> triplet : speciesTriplets) {
                    for(String taxon : triplet) {
                        out.print(taxon + " ");
                    }
                    out.println();
                }
                out.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            throw new RuntimeException("Gene trees are needed to generate triplets!");
        }
    }

    private void writeNexFiles() {
        List<String> species = new ArrayList<>(_species2alleles.keySet());
        Collections.sort(species);

        List<List<String>> triplets = new ArrayList<>();

        if(_tripletFile != null) {
            try {
                BufferedReader br = new BufferedReader(new FileReader(_tripletFile));
                String s;
                while((s = br.readLine()) != null ) {
                    s = s.trim();
                    String[] ss = s.split(" ");
                    triplets.add(Arrays.asList(ss));
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            for(int i1 = 0 ; i1 < species.size() ; i1 ++) {
                for(int i2 = i1 + 1 ; i2 < species.size() ; i2++) {
                    for(int i3 = i2 + 1 ; i3 < species.size() ; i3++) {
                        List<String> triplet = new ArrayList<>();
                        triplet.add(species.get(i1));
                        triplet.add(species.get(i2));
                        triplet.add(species.get(i3));

                        triplets.add(triplet);
                    }
                }
            }
        }

        int index = 0;
        for(List<String> triplet : triplets) {
            String cmd = _cmd;

            String tmString = "<";
            for(String s : triplet) {
                tmString += s + ":";
                for(String a : _species2alleles.get(s)) {
                    tmString += a + ",";
                }
                tmString = tmString.substring(0, tmString.length() - 1);
                tmString += ";";
            }
            tmString = tmString.substring(0, tmString.length() - 1);
            tmString += ">";

            if(cmd.contains("-tm")) {
                int tmindex = cmd.indexOf("-tm");
                int bracketBegin = cmd.indexOf("<", tmindex);
                int bracketEnd = cmd.indexOf(">", tmindex);
                cmd = cmd.substring(0, bracketBegin) + tmString + cmd.substring(bracketEnd + 1);
            }

            try {
                BufferedReader br = new BufferedReader(new FileReader(_nexTemplateFile));
                int dotPos = _nexTemplateFile.lastIndexOf(".");
                String newfilename = _nexTemplateFile.substring(0, dotPos);
                newfilename = newfilename + "_" + index + ".nex";
                PrintWriter out = new PrintWriter(newfilename);
                System.out.println("Writing " + newfilename);

                String s;
                boolean begin = false;


                while ((s = br.readLine()) != null) {


                    out.println(s);
                    if (s.replaceAll("\\s+", "").toLowerCase().startsWith("beginphylonet;")) {
                        begin = true;
                        out.println(cmd);
                        if(!cmd.trim().endsWith(";")) out.println(";");
                        out.println("END;");
                        break;
                    }
                }
                out.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            index++;
        }
    }

    @Override
    protected String produceResult() {
        System.out.println("");

        if(_mode.toLowerCase().equals("triplets")) {
            writeTriplets();
        } else if(_mode.toLowerCase().equals("nex")) {
            _cmd = parseCmd(_nexTemplateFile);
            writeNexFiles();
        } else if(_mode.equals("Result")) {
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

            if(summary == null) {
                System.out.println("Infer one more batch of trinets.");
            } else {
                Network inferred = summary.inferredNetwork;
                System.out.println("Final network:");
                System.out.println(inferred);
            }
        }





        return "";
    }
}
