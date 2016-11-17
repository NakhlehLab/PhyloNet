package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.Program;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.MC3Core;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;

import java.io.File;
import java.util.*;

/**
 * Created by dw20 on 10/26/16.
 */
@CommandName("MCMC_SEQ")
public class MCMC_SEQ extends CommandBaseFileOutMultilocusData {

    List<Alignment> alignments = new ArrayList<>();

    public MCMC_SEQ(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                    Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                    Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader, Random rand){
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams(){
        return 0;
    }

    @Override
    protected int getMaxNumParams(){
        return 40;
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        // ----- DNA sequences -----
        this.parseMultiLociData(Program.inputNexusFileName);
        if(this.sourceIdentToMultilocusData == null) {
            throw new RuntimeException("Expected data for inference! See format here: " +
                    ""); // TODO
        }
        ParamExtractor dataParam = new ParamExtractor("loci", this.params, this.errorDetected);
        if(dataParam.ContainsSwitch){
            if(dataParam.PostSwitchParam != null) {
                try {
                    if(!(dataParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    ParameterIdentList loci = (ParameterIdentList) dataParam.PostSwitchParam;
                    for(String locus: loci.Elements){
                        noError = noError && this.assertDataExists(locus, loci.getLine(), loci.getColumn());
                        if (noError) {
                            alignments.add(new Alignment(this.sourceIdentToMultilocusData.get(locus), locus));
                        }
                    }
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized data " + dataParam.PostSwitchValue,
                            dataParam.PostSwitchParam.getLine(), dataParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -loci." + dataParam.PostSwitchValue,
                        dataParam.PostSwitchParam.getLine(), dataParam.PostSwitchParam.getColumn());
            }
        } else {
            for(String key : this.sourceIdentToMultilocusData.keySet()) {
                alignments.add(new Alignment(this.sourceIdentToMultilocusData.get(key), key));
            }
        }

        // ----- MCMC Settings -----
        // chain length
        ParamExtractor clParam = new ParamExtractor("cl", this.params, this.errorDetected);
        if(clParam.ContainsSwitch){
            if(clParam.PostSwitchParam != null) {
                try {
                    Utils._CHAIN_LEN = Long.parseLong(clParam.PostSwitchValue);
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
                    Utils._BURNIN_LEN = Long.parseLong(blParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized chain length " + blParam.PostSwitchValue,
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
                    Utils._SAMPLE_FREQUENCY = Long.parseLong(sfParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized burn-in length " + sfParam.PostSwitchValue,
                            sfParam.PostSwitchParam.getLine(), sfParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -bl.",
                        sfParam.SwitchParam.getLine(), sfParam.SwitchParam.getColumn());
            }
        }

        // seed
        ParamExtractor sdParam = new ParamExtractor("sd", this.params, this.errorDetected);
        if(sdParam.ContainsSwitch){
            if(sdParam.PostSwitchParam != null) {
                try {
                    Utils._SEED = Long.parseLong(sdParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized random seed " + sdParam.PostSwitchValue,
                            sdParam.PostSwitchParam.getLine(), sdParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sd.",
                        sdParam.SwitchParam.getLine(), sdParam.SwitchParam.getColumn());
            }
        }

        // parallel threads
        ParamExtractor plParam = new ParamExtractor("pl", this.params, this.errorDetected);
        if(plParam.ContainsSwitch){
            if(plParam.PostSwitchParam != null) {
                try {
                    Utils._NUM_THREADS = Integer.parseInt(plParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized number of thread " + plParam.PostSwitchValue,
                            plParam.PostSwitchParam.getLine(), plParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -pl.",
                        plParam.SwitchParam.getLine(), plParam.SwitchParam.getColumn());
            }
        }

        // output directory
        ParamExtractor outDirParam = new ParamExtractor("outDir", this.params, this.errorDetected);
        if(outDirParam.ContainsSwitch) {
            if(outDirParam.PostSwitchParam != null) {
                try {
                    Utils._OUT_DIRECTORY = outDirParam.PostSwitchValue;
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized output directory " + outDirParam.PostSwitchValue,
                            outDirParam.PostSwitchParam.getLine(), outDirParam.PostSwitchParam.getColumn());
                }
                File outDir = new File(Utils._OUT_DIRECTORY);
                if(!outDir.exists() || !outDir.isDirectory()) {
                    errorDetected.execute("Output directory doesn't exist: " + outDirParam.PostSwitchValue,
                            outDirParam.PostSwitchParam.getLine(), outDirParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -outDir.",
                        outDirParam.SwitchParam.getLine(), outDirParam.SwitchParam.getColumn());
            }
        }

        // ----- MC3 Settings -----
        // temperatures
        ParamExtractor tpParam = new ParamExtractor("mc3", this.params, this.errorDetected);
        if(tpParam.ContainsSwitch) {
            if(tpParam.PostSwitchParam != null) {
                try {
                    if(!(tpParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    List<Double> _temperatures = new ArrayList<>();
                    ParameterIdentList temps = (ParameterIdentList) tpParam.PostSwitchParam;
                    for(String tExp: temps.Elements){
                        double t = Double.parseDouble(tExp.trim());
                        if(t < 0) {
                            throw new RuntimeException();
                        }
                        _temperatures.add(t);
                    }
                    Utils._MC3_CHAINS = _temperatures;
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -tp.",
                            tpParam.PostSwitchParam.getLine(), tpParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -tp.",
                        tpParam.SwitchParam.getLine(), tpParam.SwitchParam.getColumn());
            }
        }

        // ----- Inference Settings -----
        // maximum reticulation
        ParamExtractor mrParam = new ParamExtractor("mr", this.params, this.errorDetected);
        if(mrParam.ContainsSwitch) {
            if(mrParam.PostSwitchParam != null) {
                try {
                    Utils._NET_MAX_RETI = Integer.parseInt(mrParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized value of maximum number of reticulation nodes " + mrParam.PostSwitchValue,
                            mrParam.PostSwitchParam.getLine(), mrParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -mr.",
                        mrParam.SwitchParam.getLine(), mrParam.SwitchParam.getColumn());
            }
        }

        // taxon map
        ParamExtractor tmParam = new ParamExtractor("tm", this.params, this.errorDetected);
        if(tmParam.ContainsSwitch){
            ParamExtractorAllelListMap aaParam = new ParamExtractorAllelListMap("tm", this.params, this.errorDetected);
            noError = noError && aaParam.IsValidMap;
            if(aaParam.IsValidMap){
                Utils._TAXON_MAP = aaParam.ValueMap;
            }
        }

        // fix population size
        ParamExtractor fixPsParam = new ParamExtractor("fixps", this.params, this.errorDetected);
        if(fixPsParam.ContainsSwitch){
            if(fixPsParam.PostSwitchParam != null) {
                try {
                    Utils._ESTIMATE_POP_SIZE = false;
                    Utils._CONST_POP_SIZE = true;
                    Utils._POP_SIZE_MEAN = Double.parseDouble(fixPsParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized population size " + fixPsParam.PostSwitchValue,
                            fixPsParam.PostSwitchParam.getLine(), fixPsParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -fixps.",
                        fixPsParam.SwitchParam.getLine(), fixPsParam.SwitchParam.getColumn());
            }
        }

        // vary population sizes
        ParamExtractor varyPsParam = new ParamExtractor("varyps", this.params, this.errorDetected);
        if(varyPsParam.ContainsSwitch) {
            Utils._ESTIMATE_POP_SIZE = true;
            Utils._CONST_POP_SIZE = false;
        }

        // ----- Prior Settings -----
        // poisson parameter
        ParamExtractor ppParam = new ParamExtractor("pp", this.params, this.errorDetected);
        if(ppParam.ContainsSwitch){
            if(ppParam.PostSwitchParam != null) {
                try {
                    Utils._POISSON_PARAM = Double.parseDouble(ppParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized Poisson parameter " + ppParam.PostSwitchValue,
                            ppParam.PostSwitchParam.getLine(), ppParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -pp.",
                        ppParam.SwitchParam.getLine(), ppParam.SwitchParam.getColumn());
            }
        }

        // disable diameter prior
        ParamExtractor ddParam = new ParamExtractor("dd", this.params, this.errorDetected);
        if(ddParam.ContainsSwitch) {
            Utils._DIAMETER_PRIOR = false;
        }

        // enable time prior
        ParamExtractor eeParam = new ParamExtractor("ee", this.params, this.errorDetected);
        if(eeParam.ContainsSwitch) {
            Utils._TIMES_EXP_PRIOR = true;
        }

        // ----- Substitution Model Settings -----
        // GTR
        ParamExtractor gtrParam = new ParamExtractor("gtr", this.params, this.errorDetected);
        if(gtrParam.ContainsSwitch) {
            if(gtrParam.PostSwitchParam != null) {
                try {
                    if(!(gtrParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    double[] baseFreqs = new double[4];
                    double[] transRates = new double[6];
                    ParameterIdentList rates = (ParameterIdentList) gtrParam.PostSwitchParam;
                    int idx = 0;
                    for(String item: rates.Elements){
                        double r = Double.parseDouble(item.trim());
                        if(r < 0) {
                            throw new RuntimeException();
                        }
                        if(idx >= 4) {
                            transRates[idx - 4] = r;
                        } else {
                            baseFreqs[idx] = r;
                        }
                        idx++;
                    }
                    if(idx != 10) {
                        throw new RuntimeException();
                    }
                    Utils._BASE_FREQS = baseFreqs;
                    Utils._TRANS_RATES = transRates;
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -gtr.",
                            gtrParam.PostSwitchParam.getLine(), gtrParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -gtr.",
                        gtrParam.SwitchParam.getLine(), gtrParam.SwitchParam.getColumn());
            }
        }

        // ----- Site Model -----
        // mutation rate
        ParamExtractor muParam = new ParamExtractor("mu", this.params, this.errorDetected);
        if(muParam.ContainsSwitch){
            if(muParam.PostSwitchParam != null) {
                try {
                    Utils._MUTATION_RATE = Double.parseDouble(muParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized mutation rate " + muParam.PostSwitchValue,
                            muParam.PostSwitchParam.getLine(), muParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -mu.",
                        muParam.SwitchParam.getLine(), muParam.SwitchParam.getColumn());
            }
        }

        // ----- Starting State Settings -----
        // starting gene trees
        ParamExtractor sgtParam = new ParamExtractor("sgt", this.params, this.errorDetected);
        if(sgtParam.ContainsSwitch){
            if(sgtParam.PostSwitchParam != null) {
                try {
                    if(!(sgtParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    List<String> geneTrees = new ArrayList<>();
                    ParameterIdentList gts = (ParameterIdentList) sgtParam.PostSwitchParam;
                    for(String ident: gts.Elements){
                        noError = noError && this.assertNetworkExists(ident,
                                sgtParam.PostSwitchParam.getLine(), sgtParam.PostSwitchParam.getColumn());
                        if (noError) {
                            NetworkNonEmpty gt = this.sourceIdentToNetwork.get(ident);
                            geneTrees.add(NetworkTransformer.toENewickTree(gt));
                        }
                    }
                    Utils._START_GT_LIST = geneTrees;
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -sgt.",
                            sgtParam.PostSwitchParam.getLine(), sgtParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sgt.",
                        sgtParam.SwitchParam.getLine(), sgtParam.SwitchParam.getColumn());
            }
        }

        // starting networks
        ParamExtractor snParam = new ParamExtractor("snet", this.params, this.errorDetected);
        if(snParam.ContainsSwitch){
            if(snParam.PostSwitchParam != null) {
                try {
                    noError = noError && this.assertNetworkExists(snParam.PostSwitchValue,
                            snParam.PostSwitchParam.getLine(), snParam.PostSwitchParam.getColumn());
                    if (noError) {
                        NetworkNonEmpty net = this.sourceIdentToNetwork.get(snParam.PostSwitchValue);
                        noError = noError && net != null;
                        if (noError) {
                            Utils._START_NET = (new NetworkFactoryFromRNNetwork()).makeNetwork(net).toString();
                        }
                    }
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -snet.",
                            snParam.PostSwitchParam.getLine(), snParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -snet.",
                        snParam.SwitchParam.getLine(), snParam.SwitchParam.getColumn());
            }
        }

        // starting population size
        ParamExtractor spsParam = new ParamExtractor("sps", this.params, this.errorDetected);
        if(spsParam.ContainsSwitch){
            if(spsParam.PostSwitchParam != null) {
                try {
                    Utils._POP_SIZE_MEAN = Double.parseDouble(spsParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized population size " + spsParam.PostSwitchValue,
                            spsParam.PostSwitchParam.getLine(), spsParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sps.",
                        spsParam.SwitchParam.getLine(), spsParam.SwitchParam.getColumn());
            }
        }

        noError = noError && checkForUnknownSwitches(
                "loci",
                "cl", "bl", "sf", "sd", "pl", "outDir",
                "mc3", "mr", "tm", "fixps", "varyps",
                "pp", "dd", "ee", "gtr", "mu",
                "sgt", "snet", "sps"
        );
        checkAndSetOutFile(
                dataParam,
                clParam, blParam, sfParam, sdParam, plParam, outDirParam,
                tpParam, mrParam, tmParam, fixPsParam, varyPsParam,
                ppParam, ddParam, eeParam, gtrParam, muParam,
                sgtParam, snParam, spsParam
        );

        return  noError;
    }

    @Override
    protected String produceResult() {

        StringBuffer result = new StringBuffer("\n");
        System.out.println();

        long startTime = System.currentTimeMillis();

        MC3Core mc3 = new MC3Core(alignments);
        mc3.run();

        result.append(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));

        return result.toString();
    }

}
