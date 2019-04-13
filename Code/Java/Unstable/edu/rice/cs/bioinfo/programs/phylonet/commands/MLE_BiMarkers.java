package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.Program;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core.SimulatedAnnealing;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.Algorithms;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.RPattern;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.*;

/**
 * Created by zhujiafan on 10/26/16.
 */
@CommandName("MLE_BiMarkers")
public class MLE_BiMarkers extends CommandBaseFileOutMatrix {

    private Map<String, String> _sequence = null;

    private Double _pi0 = null;
    private boolean _diploid = false;
    private Integer _polyploid = null;
    private Character _dominant = null;
    private String _trueNetwork = null;
    private int _maxReticulations = 4;

    private int _maxExamCount = 50000;
    private int _numRuns = 100;
    private int _numOpt = 5;
    private int _maxFailure = 100;


    public MLE_BiMarkers(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
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
        return 50;
    }

    @Override
    protected boolean checkParamsForCommand(){
        boolean noError = true;

        // ------ approximated splitting -------
        ParamExtractor approxParam = new ParamExtractor("approx", this.params, this.errorDetected);
        if(approxParam.ContainsSwitch) {
            Algorithms.SWITCH_APPROX_SPLIT = true;
        }

        // ----- DNA sequences -----
        this.parseMatrixData(Program.inputNexusFileName);
        if(this.sourceIdentToMatrixData == null) {
            throw new RuntimeException("Expected data for inference! See format here: " +
                    "https://wiki.rice.edu/confluence/display/PHYLONET/");
        }

        // approximate bayesian
        ParamExtractor abcParam = new ParamExtractor("abc", this.params, this.errorDetected);
        if(abcParam.ContainsSwitch){
            SNAPPLikelihood.useApproximateBayesian = true;
        } else {
            SNAPPLikelihood.useApproximateBayesian = false;
        }

        // pseudo likelihood
        ParamExtractor pseudoParam = new ParamExtractor("pseudo", this.params, this.errorDetected);
        if(pseudoParam.ContainsSwitch){
            SNAPPLikelihood.usePseudoLikelihood = true;
        } else {
            SNAPPLikelihood.usePseudoLikelihood = false;
        }

        // ----- diploid -----
        ParamExtractor diploidParam = new ParamExtractor("diploid", this.params, this.errorDetected);
        if(diploidParam.ContainsSwitch) {
            _diploid = true;
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
                    Algorithms.HAS_DOMINANT_MARKER = true;
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized marker " + dominantParam.PostSwitchValue,
                            dominantParam.PostSwitchParam.getLine(), dominantParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -dominant." + dominantParam.PostSwitchValue,
                        dominantParam.PostSwitchParam.getLine(), dominantParam.PostSwitchParam.getColumn());
            }
        }

        // only polymorphic
        ParamExtractor onlyPolyParam = new ParamExtractor("op", this.params, this.errorDetected);
        if(onlyPolyParam.ContainsSwitch){
            SNAPPLikelihood.useOnlyPolymorphic = true;
        } else {
            SNAPPLikelihood.useOnlyPolymorphic = false;
        }

        // ----- selected taxa -----
        ParamExtractor dataParam = new ParamExtractor("taxa", this.params, this.errorDetected);
        if(dataParam.ContainsSwitch){
//            if(dataParam.PostSwitchParam != null) {
//                try {
//                    if(!(dataParam.PostSwitchParam instanceof ParameterIdentList)){
//                        throw new RuntimeException();
//                    }
//                    ParameterIdentList taxa = (ParameterIdentList) dataParam.PostSwitchParam;
//                    _sequence = new HashMap<>();
//                    for(String taxon: taxa.Elements){
//                        noError = noError && this.assertDataExists(taxon, taxa.getLine(), taxa.getColumn());
//                        if (noError) {
//                            _sequence.put(taxon, this.sourceIdentToMatrixData.get(taxon));
//                        }
//                    }
//                } catch(NumberFormatException e) {
//                    errorDetected.execute("Unrecognized data " + dataParam.PostSwitchValue,
//                            dataParam.PostSwitchParam.getLine(), dataParam.PostSwitchParam.getColumn());
//                }
//            } else {
//                errorDetected.execute("Expected value after switch -taxa." + dataParam.PostSwitchValue,
//                        dataParam.PostSwitchParam.getLine(), dataParam.PostSwitchParam.getColumn());
//            }
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

        // ----- ML Settings -----
        ParamExtractor mnoParam = new ParamExtractor("mno", this.params, this.errorDetected);
        if(mnoParam.ContainsSwitch) {
            if(mnoParam.PostSwitchParam != null) {
                try {
                    _numOpt = Integer.parseInt(mnoParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized value of maximum number of optimums " + mnoParam.PostSwitchValue,
                            mnoParam.PostSwitchParam.getLine(), mnoParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -mno.",
                        mnoParam.SwitchParam.getLine(), mnoParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor mnrParam = new ParamExtractor("mnr", this.params, this.errorDetected);
        if(mnrParam.ContainsSwitch) {
            if(mnrParam.PostSwitchParam != null) {
                try {
                    _numRuns = Integer.parseInt(mnrParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized value of maximum number of runs " + mnrParam.PostSwitchValue,
                            mnrParam.PostSwitchParam.getLine(), mnrParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -mnr.",
                        mnrParam.SwitchParam.getLine(), mnrParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor mfParam = new ParamExtractor("mf", this.params, this.errorDetected);
        if(mfParam.ContainsSwitch) {
            if(mfParam.PostSwitchParam != null) {
                try {
                    _maxFailure = Integer.parseInt(mfParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized value of maximum number of failures " + mfParam.PostSwitchValue,
                            mfParam.PostSwitchParam.getLine(), mfParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -mf.",
                        mfParam.SwitchParam.getLine(), mfParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor mecParam = new ParamExtractor("mec", this.params, this.errorDetected);
        if(mecParam.ContainsSwitch) {
            if(mecParam.PostSwitchParam != null) {
                try {
                    _maxExamCount = Integer.parseInt(mecParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized value of maximum number of examinations " + mecParam.PostSwitchValue,
                            mecParam.PostSwitchParam.getLine(), mecParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -mec.",
                        mecParam.SwitchParam.getLine(), mecParam.SwitchParam.getColumn());
            }
        }

        // ----- Inference Settings -----
        // maximum reticulation
        ParamExtractor mrParam = new ParamExtractor("mr", this.params, this.errorDetected);
        if(mrParam.ContainsSwitch) {
            if(mrParam.PostSwitchParam != null) {
                try {
                    _maxReticulations = Integer.parseInt(mrParam.PostSwitchValue);
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
                _sequence = new HashMap<>();

                Set<String> alleleSet = new HashSet<>();
                for(String species :Utils._TAXON_MAP.keySet() ) {
                    alleleSet.addAll(Utils._TAXON_MAP.get(species));
                }
                for(String alleleName : alleleSet) {
                    if(alleleSet.contains(alleleName)) {
                        _sequence.put(alleleName, this.sourceIdentToMatrixData.get(alleleName));
                    }
                }
            }
        } else {
            _sequence = new HashMap<>();
            for(String alleleName : this.sourceIdentToMatrixData.keySet()) {
                _sequence.put(alleleName, this.sourceIdentToMatrixData.get(alleleName));
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
        } else {
            Utils._CONST_POP_SIZE = true;
            Utils._ESTIMATE_POP_SIZE = true;
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

        ParamExtractor thetaWinParam = new ParamExtractor("thetawindow", this.params, this.errorDetected);
        if(thetaWinParam.ContainsSwitch){
            if(thetaWinParam.PostSwitchParam != null) {
                try {
                    Utils._POP_SIZE_WINDOW_SIZE = Double.parseDouble(thetaWinParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized theta window size " + thetaWinParam.PostSwitchValue,
                            thetaWinParam.PostSwitchParam.getLine(), thetaWinParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -thetawindow.",
                        thetaWinParam.SwitchParam.getLine(), thetaWinParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor timeWinParam = new ParamExtractor("timewindow", this.params, this.errorDetected);
        if(timeWinParam.ContainsSwitch){
            if(timeWinParam.PostSwitchParam != null) {
                try {
                    Utils._TIME_WINDOW_SIZE = Double.parseDouble(timeWinParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized time window size " + timeWinParam.PostSwitchValue,
                            timeWinParam.PostSwitchParam.getLine(), timeWinParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -timewindow.",
                        timeWinParam.SwitchParam.getLine(), timeWinParam.SwitchParam.getColumn());
            }
        }

        ParamExtractor dcParam = new ParamExtractor("dc", this.params, this.errorDetected);
        if(dcParam.ContainsSwitch){
            if(dcParam.PostSwitchParam != null) {
                try {
                    Utils.DIMENSION_CHANGE_WEIGHT = Double.parseDouble(dcParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized dimension changing rate " + dcParam.PostSwitchValue,
                            dcParam.PostSwitchParam.getLine(), dcParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -dc.",
                        dcParam.SwitchParam.getLine(), dcParam.SwitchParam.getColumn());
            }
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
                "diploid", "dominant","polyploid",
                "taxa", "op",
                "cl", "bl", "sf", "sd", "pl",
                "mc3", "mr", "tm", "fixtheta", /*"varytheta",*/
                "pp", "dd", "ee", "pi0",
                "esptheta", "snet", "truenet", "ptheta",
                "thetawindow", "timewindow", "dc",
                "pseudo", "abc","approx",
                "mno", "mnr", "mf", "mec"
        );
        checkAndSetOutFile(
                diploidParam,polyParam,
                dominantParam,
                dataParam, onlyPolyParam,
                sdParam, plParam,
                mrParam, tmParam, fixPsParam, varyPsParam,
                ppParam, ddParam, eeParam, estimatePThetaParam,
                pthetaParam, snParam, tnParam, pi0Param,
                thetaWinParam, timeWinParam, dcParam,
                pseudoParam, abcParam,approxParam,
                mnoParam, mnrParam, mfParam, mecParam
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
                    if(s.charAt(i) != '?') {
                        count[s.charAt(i) - '0']++;
                    }
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

        System.out.println("PI0 = " + pi[0] + " PI1 = " + pi[1]);

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

        Utils.printSettings();

        List<MarkerSeq> alnwarp = new ArrayList<>();
        alnwarp.add(new MarkerSeq(_sequence));

        if(_polyploid != null) {
            alnwarp.get(0)._RPatterns = SNAPPLikelihood.polyploidSequenceToPatterns(allele2species, alnwarp, _polyploid);
        } else if(_diploid && _dominant == null) {
            alnwarp.get(0)._RPatterns = SNAPPLikelihood.diploidSequenceToPatterns(allele2species, alnwarp);
        } else {
            alnwarp.get(0)._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(allele2species, alnwarp);
        }

        // TODO: fix (pseudo likelihood) for polyploid
        alnwarp.get(0)._diploid = _diploid;
        alnwarp.get(0)._dominant = _dominant;
        alnwarp.get(0)._polyploid = _polyploid;

        double polyCount = 0.0;
        for(RPattern key : alnwarp.get(0)._RPatterns.keySet()) {
            if(!key.isMonomorphic()) {
                polyCount += alnwarp.get(0)._RPatterns.get(key)[0];
            }
        }
        System.out.println("Polymorphic sites: " + polyCount);

        if(_dominant != null) {
            Algorithms.HAS_DOMINANT_MARKER = true;
        } else {
            Algorithms.HAS_DOMINANT_MARKER = false;
        }

        if(_trueNetwork != null) {
            System.out.println("True Network: " + _trueNetwork);
            double trueRootPopSize = Double.parseDouble(_trueNetwork.substring(_trueNetwork.indexOf('[') + 1, _trueNetwork.indexOf(']')));
            Network cloneNetwork = edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(_trueNetwork.substring(_trueNetwork.indexOf(']') + 1));
            cloneNetwork.getRoot().setRootPopSize(trueRootPopSize);
            if(!SNAPPLikelihood.usePseudoLikelihood && !SNAPPLikelihood.useApproximateBayesian){
                System.out.println("True Likelihood = " + SNAPPLikelihood.computeSNAPPLikelihoodST(cloneNetwork, alnwarp.get(0)._RPatterns, BAGTRModel));
            }

            if(SNAPPLikelihood.usePseudoLikelihood) {
                System.out.println("True Pseudo-Likelihood = " + SNAPPLikelihood.computeSNAPPPseudoLikelihood(cloneNetwork, allele2species, alnwarp, BAGTRModel));
            }

            if(SNAPPLikelihood.useApproximateBayesian) {
                System.out.println("True Approximate Bayesian = " + SNAPPLikelihood.computeApproximateBayesian(cloneNetwork, allele2species, alnwarp, BAGTRModel, new HashMap<>()));
            }
        }

        Utils._NET_MAX_RETI = _maxReticulations;
        Utils._MCMC = false;
        SimulatedAnnealing sa = new SimulatedAnnealing();
        sa.setSeed(Utils._SEED);
        LinkedList<Tuple<Network,Double>> resultList = new LinkedList<>();
        sa.search(alnwarp, BAGTRModel, _numOpt, _numRuns, _maxExamCount, _maxFailure, false, resultList);

        result.append(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));

        return result.toString();
    }

}
