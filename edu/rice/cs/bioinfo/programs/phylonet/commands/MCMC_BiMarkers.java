package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.Program;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core.MC3Core;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;

import java.io.File;
import java.util.*;

/**
 * Created by zhujiafan on 10/26/16.
 */
@CommandName("MCMC_BiMarkers")
public class MCMC_BiMarkers extends CommandBaseFileOutMatrix {

    private Map<String, String> _sequence = null;

    private Double _pi0 = null;
    private boolean _diploid = false;
    private Character _dominant = null;

    public MCMC_BiMarkers(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
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
        this.parseMatrixData(Program.inputNexusFileName);
        if(this.sourceIdentToMatrixData == null) {
            throw new RuntimeException("Expected data for inference! See format here: " +
                    "https://wiki.rice.edu/confluence/display/PHYLONET/");
        }

        // ----- diploid -----
        ParamExtractor diploidParam = new ParamExtractor("diploid", this.params, this.errorDetected);
        if(diploidParam.ContainsSwitch) {
            _diploid = true;
        }

        ParamExtractor dominantParam = new ParamExtractor("dominant", this.params, this.errorDetected);
        if(dominantParam.ContainsSwitch) {
            if(!_diploid) {
                errorDetected.execute("Only diploid has dominant markers.",
                        dominantParam.PostSwitchParam.getLine(), dominantParam.PostSwitchParam.getColumn());
            }

            if(dominantParam.PostSwitchParam != null) {
                try {
                    int marker = Integer.parseInt(dominantParam.PostSwitchValue);
                    if(marker != 0 && marker != 1)
                        throw new NumberFormatException();
                    _dominant = (char)('0' + marker);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized marker " + dominantParam.PostSwitchValue,
                            dominantParam.PostSwitchParam.getLine(), dominantParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -dominant." + dominantParam.PostSwitchValue,
                        dominantParam.PostSwitchParam.getLine(), dominantParam.PostSwitchParam.getColumn());
            }
        }

        // ----- selected taxa -----
        ParamExtractor dataParam = new ParamExtractor("taxa", this.params, this.errorDetected);
        if(dataParam.ContainsSwitch){
            if(dataParam.PostSwitchParam != null) {
                try {
                    if(!(dataParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    ParameterIdentList taxa = (ParameterIdentList) dataParam.PostSwitchParam;
                    _sequence = new HashMap<>();
                    for(String taxon: taxa.Elements){
                        noError = noError && this.assertDataExists(taxon, taxa.getLine(), taxa.getColumn());
                        if (noError) {
                            _sequence.put(taxon, this.sourceIdentToMatrixData.get(taxon));
                        }
                    }
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized data " + dataParam.PostSwitchValue,
                            dataParam.PostSwitchParam.getLine(), dataParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -taxa." + dataParam.PostSwitchValue,
                        dataParam.PostSwitchParam.getLine(), dataParam.PostSwitchParam.getColumn());
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

        // fix theta
        ParamExtractor fixPsParam = new ParamExtractor("fixtheta", this.params, this.errorDetected);
        if(fixPsParam.ContainsSwitch){
            if(fixPsParam.PostSwitchParam != null) {
                try {
                    Utils._ESTIMATE_POP_SIZE = false;
                    Utils._CONST_POP_SIZE = true;
                    Utils._POP_SIZE_MEAN = Double.parseDouble(fixPsParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized theta " + fixPsParam.PostSwitchValue,
                            fixPsParam.PostSwitchParam.getLine(), fixPsParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -fixtheta.",
                        fixPsParam.SwitchParam.getLine(), fixPsParam.SwitchParam.getColumn());
            }
        }

        // varying theta in different branch
        ParamExtractor varyPsParam = new ParamExtractor("varytheta", this.params, this.errorDetected);
        if(varyPsParam.ContainsSwitch) {
            Utils._ESTIMATE_POP_SIZE = true;
            Utils._CONST_POP_SIZE = false;
        }

        // estimate theta prior
        ParamExtractor estimatePThetaParam = new ParamExtractor("esptheta", this.params, this.errorDetected);
        if(estimatePThetaParam.ContainsSwitch) {
            Utils._ESTIMATE_POP_PARAM = true;
        } else {
            Utils._ESTIMATE_POP_PARAM = false;
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

        // enable diameter prior
        ParamExtractor ddParam = new ParamExtractor("dd", this.params, this.errorDetected);
        if(ddParam.ContainsSwitch) {
            Utils._DIAMETER_PRIOR = true;
        } else {
            Utils._DIAMETER_PRIOR = false;
        }

        // enable time prior
        ParamExtractor eeParam = new ParamExtractor("ee", this.params, this.errorDetected);
        if(eeParam.ContainsSwitch) {
            Utils._TIMES_EXP_PRIOR = true;
            if(eeParam.PostSwitchParam != null) {
                try {
                    Utils.EXP_PARAM = Double.parseDouble(eeParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized exp parameter " + eeParam.PostSwitchValue,
                            eeParam.PostSwitchParam.getLine(), eeParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -ee.",
                        eeParam.SwitchParam.getLine(), eeParam.SwitchParam.getColumn());
            }
        } else {
            Utils._TIMES_EXP_PRIOR = false;
        }

        // ----- Starting State Settings -----

        // starting network
        ParamExtractor snParam = new ParamExtractor("snet", this.params, this.errorDetected);
        if(snParam.ContainsSwitch){
            if(snParam.PostSwitchParam != null) {
                if (noError) {
                    Utils._START_NET = snParam.PostSwitchValue;
                }
            } else {
                errorDetected.execute("Expected value after switch -snet.",
                        snParam.SwitchParam.getLine(), snParam.SwitchParam.getColumn());
            }
        }

        // starting population size
        ParamExtractor pthetaParam = new ParamExtractor("ptheta", this.params, this.errorDetected);
        if(pthetaParam.ContainsSwitch){
            if(pthetaParam.PostSwitchParam != null) {
                try {
                    Utils._POP_SIZE_MEAN = Double.parseDouble(pthetaParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized theta prior " + pthetaParam.PostSwitchValue,
                            pthetaParam.PostSwitchParam.getLine(), pthetaParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -ptheta.",
                        pthetaParam.SwitchParam.getLine(), pthetaParam.SwitchParam.getColumn());
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

        noError = noError && checkForUnknownSwitches(
                "diploid", "dominant",
                "taxa",
                "cl", "bl", "sf", "sd", "pl",
                "mc3", "mr", "tm", "fixtheta", "varytheta",
                "pp", "dd", "ee", "pi0",
                "esptheta", "snet", "ptheta"
        );
        checkAndSetOutFile(
                diploidParam,
                dominantParam,
                dataParam,
                clParam, blParam, sfParam, sdParam, plParam,
                tpParam, mrParam, tmParam, fixPsParam, varyPsParam,
                ppParam, ddParam, eeParam, estimatePThetaParam,
                pthetaParam, snParam, pi0Param
        );

        return  noError;
    }

    @Override
    protected String produceResult() {

        StringBuffer result = new StringBuffer("\n");
        System.out.println("\nOutput files under " + Utils._OUT_DIRECTORY);

        long startTime = System.currentTimeMillis();

        if(_dominant != null && _dominant.equals('1')) {
            if(_pi0 != null)
                _pi0 = 1.0 - _pi0;

            for (String taxon : _sequence.keySet()) {
                StringBuilder newseq = new StringBuilder(_sequence.get(taxon));
                for (int i = 0; i < newseq.length(); i++) {
                    if (newseq.charAt(i) == '0')
                        newseq.setCharAt(i, '1');
                    else if (newseq.charAt(i) == '1')
                        newseq.setCharAt(i, '0');
                }
                _sequence.put(taxon, newseq.toString());
                System.out.println(taxon + " " + _sequence.get(taxon).substring(0, 10) + "..." + _sequence.get(taxon).substring(_sequence.get(taxon).length() - 10, _sequence.get(taxon).length()));
            }
        }

        double [] pi =  new double[2];
        double [] rate = new double[1];

        if(_pi0 == null) {
            int totalSites = 0;
            int [] count = new int[3];
            count[0] = count[1] = count[2] = 0;
            for(String taxon : _sequence.keySet()) {
                String s = _sequence.get(taxon);
                for(int i = 0 ; i < s.length() ; i++) {
                    count[s.charAt(i) - '0']++;
                    totalSites++;
                }
            }

            if(_diploid && _dominant == null)
                pi[1] = 1.0 * (count[1] + 2 * count[2]) / (totalSites * 2);
            else
                pi[1] = 1.0 * count[1] / totalSites;

            pi[0] = 1 - pi[1];
            rate[0] = 1 / (2*pi[0]);
        } else {
            pi[0] = _pi0;
            pi[1] = 1 - pi[0];
            rate[0] = 1 / (2*pi[0]);
        }

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(pi, rate);

        Map<String, String> allele2species = null;
        if(Utils._TAXON_MAP != null) {
            allele2species = new HashMap<>();
            for (String species : Utils._TAXON_MAP.keySet()) {
                for (String allele : Utils._TAXON_MAP.get(species)) {
                    allele2species.put(allele, species);
                }
            }
        }

        List<Alignment> alnwarp = new ArrayList<>();
        alnwarp.add(new Alignment(_sequence));

        if(_diploid && _dominant == null) {
            alnwarp.get(0)._RPatterns = SNAPPLikelihood.diploidSequenceToPatterns(allele2species, alnwarp);
        } else {
            alnwarp.get(0)._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(allele2species, alnwarp);
        }
        MC3Core mc3 = new MC3Core(alnwarp, BAGTRModel);
        mc3.run();

        result.append(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));

        return result.toString();
    }

}
