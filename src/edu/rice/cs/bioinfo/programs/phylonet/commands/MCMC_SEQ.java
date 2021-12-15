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
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary.Convergence;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.File;
import java.util.*;

/**
 * Created by dw20 on 10/26/16.
 */
@CommandName("MCMC_SEQ")
public class MCMC_SEQ extends CommandBaseFileOutMultilocusData {

    private List<Alignment> alignments = new ArrayList<>();
    private String gtoutgroup = null;

    // for summarization
    private List<String> _files = new ArrayList<>();

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

        // summary
        ParamExtractor sumParam = new ParamExtractor("sum", this.params, this.errorDetected);
        if(sumParam.ContainsSwitch) {

            if (sumParam.PostSwitchParam != null) {
                try {
                    if (!(sumParam.PostSwitchParam instanceof ParameterIdentList)) {
                        throw new RuntimeException();
                    }
                    ParameterIdentList temps = (ParameterIdentList) sumParam.PostSwitchParam;
                    for (String tExp : temps.Elements) {
                        System.out.println(tExp);
                        _files.add(tExp);
                    }
                } catch (NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -sum.",
                            sumParam.PostSwitchParam.getLine(), sumParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sum.",
                        sumParam.SwitchParam.getLine(), sumParam.SwitchParam.getColumn());
            }

            noError = noError && checkForUnknownSwitches("sum");
            checkAndSetOutFile(sumParam);
            return noError;
        }

        // ----- DNA sequences -----
        this.parseMultiLociData(Program.inputNexusFileName);
        if(this.sourceIdentToMultilocusData == null) {
            throw new RuntimeException("Expected data for inference! See format here: " +
                    "https://wiki.rice.edu/confluence/display/PHYLONET/Multilocus+Data");
        }

        // ----- diploid phasing -----
        ParamExtractor diploidParam = new ParamExtractor("diploid", this.params, this.errorDetected);
        if(diploidParam.ContainsSwitch) {
            if(diploidParam.PostSwitchParam != null) {
                try {
                    if(!(diploidParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    Set<String> diploidSpecies = new HashSet<>();
                    ParameterIdentList taxa = (ParameterIdentList) diploidParam.PostSwitchParam;
                    for(String taxon: taxa.Elements){
                        diploidSpecies.add(taxon);
                    }
                    Utils._DIPLOID_SPECIES = diploidSpecies;
                    Utils._PHASING = diploidSpecies.size() > 0;
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -diploid.",
                            diploidParam.PostSwitchParam.getLine(), diploidParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -diploid.",
                        diploidParam.SwitchParam.getLine(), diploidParam.SwitchParam.getColumn());
            }
        }


        // taxon map
        ParamExtractor tmParam = new ParamExtractor("tm", this.params, this.errorDetected);
        if(tmParam.ContainsSwitch){
            ParamExtractorAllelListMap aaParam = new ParamExtractorAllelListMap("tm", this.params, this.errorDetected);
            noError = noError && aaParam.IsValidMap;
            if(aaParam.IsValidMap){
                Utils._TAXON_MAP = aaParam.ValueMap;

                if(Utils._TAXON_MAP != null) {
                    Set<String> alleleSet = new HashSet<>();
                    for(String species :Utils._TAXON_MAP.keySet() ) {
                        alleleSet.addAll(Utils._TAXON_MAP.get(species));
                    }
                    for(String locusName : this.sourceIdentToMultilocusData.keySet()) {
                        Map<String, String> locus = this.sourceIdentToMultilocusData.get(locusName);
                        Map<String, String> newLocus = new HashMap<>();
                        for(String alleleName : locus.keySet()) {
                            if(alleleSet.contains(alleleName)) {
                                newLocus.put(alleleName, locus.get(alleleName));
                            }
                        }
                        this.sourceIdentToMultilocusData.put(locusName, newLocus);
                    }
                }

            }
        }



        // ----- selected loci -----
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
                    Randomizer.setSeed(Utils._SEED);
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
        ParamExtractor dirParam = new ParamExtractor("dir", this.params, this.errorDetected);
        if(dirParam.ContainsSwitch) {
            if(dirParam.PostSwitchParam != null) {
                try {
                    Utils._OUT_DIRECTORY = dirParam.PostSwitchValue;
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized output directory " + dirParam.PostSwitchValue,
                            dirParam.PostSwitchParam.getLine(), dirParam.PostSwitchParam.getColumn());
                }
                File dir = new File(Utils._OUT_DIRECTORY);
                if(!dir.exists() || !dir.isDirectory()) {
                    errorDetected.execute("Output directory doesn't exist: " + dirParam.PostSwitchValue,
                            dirParam.PostSwitchParam.getLine(), dirParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -dir.",
                        dirParam.SwitchParam.getLine(), dirParam.SwitchParam.getColumn());
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

        // sample embeddings
        ParamExtractor seParam = new ParamExtractor("se", this.params, this.errorDetected);
        if(seParam.ContainsSwitch) {
            Utils.SAMPLE_EMBEDDINGS = true;
        }

        // use pseudo likelihood (Yun 2015) to compute net likelihood
        ParamExtractor pseudoParam = new ParamExtractor("pseudo", this.params, this.errorDetected);
        if(pseudoParam.ContainsSwitch) {
            Utils.PSEUDO_LIKELIHOOD = true;
        }

        // fix all gene tree topologies
        ParamExtractor gtburninParam = new ParamExtractor("gtburnin", this.params, this.errorDetected);
        if(gtburninParam.ContainsSwitch) {
            Utils._START_GT_BURN_IN = true;
        }

        // fix all gene tree topologies
        ParamExtractor fixgttopoParam = new ParamExtractor("fixgttopo", this.params, this.errorDetected);
        if(fixgttopoParam.ContainsSwitch) {
            Utils._FIX_GENE_TREE_TOPOLOGIES = true;
        }

        // fix all gene trees
        ParamExtractor fixgtsParam = new ParamExtractor("fixgts", this.params, this.errorDetected);
        if(fixgtsParam.ContainsSwitch) {
            Utils._FIX_GENE_TREES = true;
        }

        // fix net topology
        ParamExtractor fixNetTopoParam = new ParamExtractor("fixnettopo", this.params, this.errorDetected);
        if(fixNetTopoParam.ContainsSwitch) {
            Utils._FIX_NET_TOPOLOGY = true;
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

        // vary mutation rate of loci
        ParamExtractor murateParam = new ParamExtractor("murate", this.params, this.errorDetected);
        if(murateParam.ContainsSwitch) {
            Utils.SAMPLE_MUTATION_RATE = true;
        }

        //mutation rate parameter input of loci
        ParamExtractor murateListParam = new ParamExtractor("mupi", this.params, this.errorDetected);
        if(murateListParam.ContainsSwitch) {
            if(murateListParam.PostSwitchParam != null) {
                try {
                    if(!(murateListParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    ParameterIdentList rates = (ParameterIdentList) murateListParam.PostSwitchParam;
                    for(String item: rates.Elements){
                        double r = Double.parseDouble(item.trim());
                        Utils._PARAMETER_INPUT.add(r);
                    }
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -mupi.",
                            murateListParam.PostSwitchParam.getLine(), murateListParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -mupi.",
                        murateListParam.SwitchParam.getLine(), murateListParam.SwitchParam.getColumn());
            }
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
        ParamExtractor disableAllPriorParam = new ParamExtractor("disableallprior", this.params, this.errorDetected);
        if(disableAllPriorParam.ContainsSwitch) {
            Utils._DISABLE_ALL_PRIOR = true;
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
                    Map<String, String> geneTrees = new HashMap<>();
                    ParameterIdentList gts = (ParameterIdentList) sgtParam.PostSwitchParam;
                    for(String ident: gts.Elements){
                        noError = noError && this.assertNetworkExists(ident,
                                sgtParam.PostSwitchParam.getLine(), sgtParam.PostSwitchParam.getColumn());
                        if (noError) {
                            NetworkNonEmpty gt = this.sourceIdentToNetwork.get(ident);
                            geneTrees.put(ident, NetworkTransformer.toENewickTree(gt));
                        }
                    }
                    List<Map<String, String>> startGTs = new ArrayList<>();
                    startGTs.add(geneTrees);
                    Utils._START_GT_LIST = startGTs;
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -sgt.",
                            sgtParam.PostSwitchParam.getLine(), sgtParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -sgt.",
                        sgtParam.SwitchParam.getLine(), sgtParam.SwitchParam.getColumn());
            }
        }

        // If gene trees are unrooted, this can root them.
        ParamExtractor gtoutgroupParam = new ParamExtractor("gtoutgroup", this.params, this.errorDetected);
        if(gtoutgroupParam.ContainsSwitch){
            if(gtoutgroupParam.PostSwitchParam != null) {
                gtoutgroup = gtoutgroupParam.PostSwitchValue;
            } else {
                errorDetected.execute("Expected value after switch -gtoutgroup.",
                        gtoutgroupParam.SwitchParam.getLine(), gtoutgroupParam.SwitchParam.getColumn());
            }
        }

        // starting networks
        ParamExtractor snParam = new ParamExtractor("snet", this.params, this.errorDetected);
        if(snParam.ContainsSwitch){
            if(snParam.PostSwitchParam != null) {
                try {
                    if(!(snParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    List<String> snets = new ArrayList<>();
                    ParameterIdentList nets = (ParameterIdentList) snParam.PostSwitchParam;
                    NetworkFactoryFromRNNetwork factory = new NetworkFactoryFromRNNetwork();
                    for(String ident: nets.Elements){
                        noError = noError && this.assertNetworkExists(ident,
                                snParam.PostSwitchParam.getLine(), snParam.PostSwitchParam.getColumn());
                        if (noError) {
                            NetworkNonEmpty net = this.sourceIdentToNetwork.get(ident);
                            snets.add(factory.makeNetwork(net).toString());
                        }
                    }
                    Utils._START_NET = snets;
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

        // pre-burn-in
        ParamExtractor preParam = new ParamExtractor("pre", this.params, this.errorDetected);
        if(preParam.ContainsSwitch){
            if(preParam.PostSwitchParam != null) {
                try {
                    int pre = Integer.parseInt(preParam.PostSwitchValue);
                    if(pre < 0) {
                        throw new RuntimeException("pre-burn-in iterations should be an integer >= 0");
                    }
                    Utils._PRE_BURN_IN_ITER = pre;
                    Utils._PRE_BURN_IN = true;
                } catch(NumberFormatException e) {
                    errorDetected.execute(
                            "Unrecognized pre-burn-in iterations (an integer >= 0)" + preParam.PostSwitchValue,
                            preParam.PostSwitchParam.getLine(), preParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -pre.",
                        preParam.SwitchParam.getLine(), preParam.SwitchParam.getColumn());
            }
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
                    Utils._SUBSTITUTION_MODEL = "GTR";
                } catch(NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -gtr.",
                            gtrParam.PostSwitchParam.getLine(), gtrParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -gtr.",
                        gtrParam.SwitchParam.getLine(), gtrParam.SwitchParam.getColumn());
            }
        }

        noError = noError && checkForUnknownSwitches(
                "diploid",
                "loci",
                "cl", "bl", "sf", "sd", "pl", "dir",
                "mc3", "mr", "tm", "fixps", "varyps",
                "fixnettopo", "fixgts", "disableallprior",
                "fixgttopo", "gtburnin",
                "pp", "dd", "ee", "mu", "se",
                "sgt", "snet", "sps", "pre", "gtr",
                "gtoutgroup", "pseudo", "murate", "mupi"
        );
        checkAndSetOutFile(
                diploidParam,
                dataParam,
                clParam, blParam, sfParam, sdParam, plParam, dirParam,
                tpParam, mrParam, tmParam, fixPsParam, varyPsParam,
                fixNetTopoParam, fixgtsParam, disableAllPriorParam,
                fixgttopoParam, gtburninParam,
                ppParam, ddParam, eeParam, muParam, seParam,
                sgtParam, snParam, spsParam, gtrParam,
                gtoutgroupParam, pseudoParam, murateParam
        );

        return  noError;
    }

    private void preprocessGeneTrees() {
        // Remove the outgroup of gene trees from data.
        for(String species : Utils._TAXON_MAP.keySet()) {
            if(Utils._TAXON_MAP.get(species).contains(gtoutgroup)) {
                Utils._TAXON_MAP.get(species).remove(gtoutgroup);
                for (String locus : this.sourceIdentToMultilocusData.keySet()) {
                    this.sourceIdentToMultilocusData.get(locus).remove(gtoutgroup);
                }
            }
        }

        for(int i = 0 ; i < Utils._START_GT_LIST.size() ; i++) {
            for (String locus : Utils._START_GT_LIST.get(i).keySet()) {
                STITree stitree = new STITree(Trees.readTree(Utils._START_GT_LIST.get(i).get(locus)));
                List<String> leaves = new ArrayList<>();
                for(String species : Utils._TAXON_MAP.keySet()) {
                    for(String allele : Utils._TAXON_MAP.get(species)) {
                        leaves.add(allele);
                    }
                }

                stitree.rerootTreeAtEdge(gtoutgroup);
                stitree.removeNode(gtoutgroup);
                Trees.removeBinaryNodes(stitree);
                stitree.constrainByLeaves(leaves);

                for(Object nodeObj : stitree.postTraverse()) {
                    STINode node = (STINode) nodeObj;
                    node.setParentDistance(TNode.NO_DISTANCE);
                }

                Utils._START_GT_LIST.get(i).put(locus, stitree.toNewick());
            }
        }
    }

    @Override
    protected String produceResult() {

        if(_files.size() > 0) {
            System.out.println();
            Convergence conv = new Convergence(_files);
            conv.summarizeTopo(true);
//            conv.computePSRF();
//            conv.generateSRQ();
//            conv.generateTracePlot();
            return "";
        }

        if(Utils._PHASING) {
            Utils.taxonMapPhasing(this.sourceIdentToMultilocusData.values().iterator().next().keySet());
            Utils._SUBSTITUTION_MODEL = "JC";
        }

        System.out.println("\nOutput files under " + Utils._OUT_DIRECTORY);

        long startTime = System.currentTimeMillis();

        if(Utils._START_GT_LIST != null && gtoutgroup != null) {
            preprocessGeneTrees();
        }

        Collections.sort(alignments);
        MC3Core mc3 = new MC3Core(alignments);
        mc3.run();

        System.out.printf("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0);

        return "";
    }

}
