package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.Program;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.RecombinationRate;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.VariationalModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.Prior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.VariationalVariable;

import java.util.*;


/**
 * Created by Xinhao Liu on 8/25/20.
 */
@CommandName("VI_coalHMM")
public class VI_coalHMM extends CommandBaseFileOutMatrix {
    private static final String _website = "https://wiki.rice.edu/confluence/display/PHYLONET/List+of+PhyloNet+Commands";

    private String _startTree = null;
    private double _originalMu;
    private double _recombRate;

    public VI_coalHMM(SyntaxCommand motivatingCommand, ArrayList<Parameter> params,
                      Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                      Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader, Random rand) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 0;
    }

    @Override
    protected int getMaxNumParams() {
        return 40;
    }

    @Override
    protected boolean checkParamsForCommand() {
        boolean noError = true;

        // ----- DNA sequences -----
        this.parseMatrixData(Program.inputNexusFileName);
        if(this.sourceIdentToMatrixData == null) {
            throw new RuntimeException("Required data for inference! See format here: " +
                    _website);
        }
        Utils.DATA = new Alignment(this.sourceIdentToMatrixData);

        // ----- Starting State Settings -----
        // starting tree with node heights
        ParamExtractor stParam = new ParamExtractor("st", this.params, this.errorDetected);
        if (stParam.ContainsSwitch) {
            if (stParam.PostSwitchParam != null) {
                try {
                    if(!(stParam.PostSwitchParam instanceof ParameterIdentList)){
                        throw new RuntimeException();
                    }
                    ParameterIdentList trees = (ParameterIdentList) stParam.PostSwitchParam;
                    for (String ident: trees.Elements) {
                        noError = noError && this.assertNetworkExists(ident,
                                stParam.PostSwitchParam.getLine(), stParam.PostSwitchParam.getColumn());
                        if (noError) {
                            NetworkNonEmpty gt = this.sourceIdentToNetwork.get(ident);
                            _startTree = NetworkTransformer.toENewickTree(gt);
                        }
                    }
                } catch (NumberFormatException e) {
                    errorDetected.execute("Invalid value after switch -st.",
                            stParam.PostSwitchParam.getLine(), stParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected newick tree string after switch -st.",
                        stParam.SwitchParam.getLine(), stParam.SwitchParam.getColumn());
            }
        } else {
            throw new RuntimeException("Required starting tree -st for inference! See format here: " + _website);
        }

        // mutation rate in unit of expected number of mutations per site per generation
        ParamExtractor muParam = new ParamExtractor("mu", this.params, this.errorDetected);
        if (muParam.ContainsSwitch) {
            if (muParam.PostSwitchParam != null) {
                try {
                    _originalMu = Double.parseDouble(muParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized mutation rate " + muParam.PostSwitchValue,
                            muParam.PostSwitchParam.getLine(), muParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -mu.",
                        muParam.SwitchParam.getLine(), muParam.SwitchParam.getColumn());
            }
        } else {
            throw new RuntimeException("Required mutation rate -mu for inference! See format here: " +
                    _website);
        }

        // recombination rate in unit of expected number of recombinations per site per generation
        ParamExtractor rhoParam = new ParamExtractor("rho", this.params, this.errorDetected);
        if (rhoParam.ContainsSwitch) {
            if (rhoParam.PostSwitchParam != null) {
                try {
                    _recombRate = Double.parseDouble(rhoParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized recombination rate " + rhoParam.PostSwitchValue,
                            rhoParam.PostSwitchParam.getLine(), rhoParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -rho.",
                        rhoParam.SwitchParam.getLine(), rhoParam.SwitchParam.getColumn());
            }
        } else {
            throw new RuntimeException("Required recombination rate -rho for inference! See format here: " +
                    _website);
        }

        // node height variational posterior initial standard deviation, default 20000
        ParamExtractor nhsigmaParam = new ParamExtractor("nhsigma", this.params, this.errorDetected);
        if (nhsigmaParam.ContainsSwitch) {
            if (nhsigmaParam.PostSwitchParam != null) {
                try {
                    Utils.NODE_HEIGHT_INIT_STDDEV = Double.parseDouble(nhsigmaParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized node height posterior initial standard deviation " + nhsigmaParam.PostSwitchValue,
                            nhsigmaParam.PostSwitchParam.getLine(), nhsigmaParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -nhsigma.",
                        nhsigmaParam.SwitchParam.getLine(), nhsigmaParam.SwitchParam.getColumn());
            }
        }

        // population size variational posterior initial standard deviation, default 10000
        ParamExtractor pssigmaParam = new ParamExtractor("pssigma", this.params, this.errorDetected);
        if (pssigmaParam.ContainsSwitch) {
            if (pssigmaParam.PostSwitchParam != null) {
                try {
                    Utils.POP_SIZE_INIT_STDDEV = Double.parseDouble(pssigmaParam.PostSwitchValue);
                } catch(NumberFormatException e) {
                    errorDetected.execute("Unrecognized population size posterior initial standard deviation " + pssigmaParam.PostSwitchValue,
                            pssigmaParam.PostSwitchParam.getLine(), pssigmaParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -pssigma.",
                        pssigmaParam.SwitchParam.getLine(), pssigmaParam.SwitchParam.getColumn());
            }
        }

        // ----- Prior Settings -----
        // mean of the gamma prior on population sizes of each branch, default 50000
        ParamExtractor pspParam = new ParamExtractor("psp", this.params, this.errorDetected);
        if (pspParam.ContainsSwitch) {
            if (pspParam.PostSwitchParam != null) {
                try {
                    Utils.POP_SIZE_MEAN = Integer.parseInt(pspParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized population size gamma prior mean " + pspParam.PostSwitchValue,
                            pspParam.PostSwitchParam.getLine(), pspParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -psp.",
                        pspParam.SwitchParam.getLine(), pspParam.SwitchParam.getColumn());
            }
        }

        // ----- Likelihood Simulator Settings -----
        // N0 for ms, default 10000
        ParamExtractor n0Param = new ParamExtractor("n0", this.params, this.errorDetected);
        if (n0Param.ContainsSwitch) {
            if (n0Param.PostSwitchParam != null) {
                try {
                    Utils.N0 = Integer.parseInt(n0Param.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized N0 " + n0Param.PostSwitchValue,
                            n0Param.PostSwitchParam.getLine(), n0Param.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -n0.",
                        n0Param.SwitchParam.getLine(), n0Param.SwitchParam.getColumn());
            }
        }

        // -r for ms
        ParamExtractor rParam = new ParamExtractor("r", this.params, this.errorDetected);
        if (rParam.ContainsSwitch) {
            if (rParam.PostSwitchParam != null) {
                try {
                    Utils.CROSS_OVER_RATE = Integer.parseInt(rParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized -r " + rParam.PostSwitchValue,
                            rParam.PostSwitchParam.getLine(), rParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -r.",
                        rParam.SwitchParam.getLine(), rParam.SwitchParam.getColumn());
            }
        } else {
            throw new RuntimeException("Required crossover rate -r for inference! See format here: " +
                    _website);
        }

        // number of subbranches per branch
        ParamExtractor nbParam = new ParamExtractor("nb", this.params, this.errorDetected);
        if (nbParam.ContainsSwitch) {
            if (nbParam.PostSwitchParam != null) {
                try {
                    Utils.NUM_BIN = Integer.parseInt(nbParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized number of sub-branches " + nbParam.PostSwitchValue,
                            nbParam.PostSwitchParam.getLine(), nbParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -nb.",
                        nbParam.SwitchParam.getLine(), nbParam.SwitchParam.getColumn());
            }
        } else {
            throw new RuntimeException("Required num_bin -nb for inference! See format here: " +
                    _website);
        }

        // ----- BBVI Settings -----
        // number of samples per iteration, default 50
        ParamExtractor nsParam = new ParamExtractor("ns", this.params, this.errorDetected);
        if (nsParam.ContainsSwitch) {
            if (nsParam.PostSwitchParam != null) {
                try {
                    Utils.nSamples = Integer.parseInt(nsParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized number of samples per iteration " + nsParam.PostSwitchValue,
                            nsParam.PostSwitchParam.getLine(), nsParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -ns.",
                        nsParam.SwitchParam.getLine(), nsParam.SwitchParam.getColumn());
            }
        }

        // number of iterations
        ParamExtractor niterParam = new ParamExtractor("niter", this.params, this.errorDetected);
        if (niterParam.ContainsSwitch) {
            if (niterParam.PostSwitchParam != null) {
                try {
                    Utils.nIterations = Integer.parseInt(niterParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized number of iterations " + niterParam.PostSwitchValue,
                            niterParam.PostSwitchParam.getLine(), niterParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -niter.",
                        niterParam.SwitchParam.getLine(), niterParam.SwitchParam.getColumn());
            }
        }

        // learning rate for mean of node height variational posterior, default 20000
        ParamExtractor nhmeanlrParam = new ParamExtractor("nhmeanlr", this.params, this.errorDetected);
        if (nhmeanlrParam.ContainsSwitch) {
            if (nhmeanlrParam.PostSwitchParam != null) {
                try {
                    Utils.NODE_HEIGHT_MEAN_LEARNING_RATE = Double.parseDouble(nhmeanlrParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized learning rate for mean of node height posterior " + nhmeanlrParam.PostSwitchValue,
                            nhmeanlrParam.PostSwitchParam.getLine(), nhmeanlrParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -nhmeanlr.",
                        nhmeanlrParam.SwitchParam.getLine(), nhmeanlrParam.SwitchParam.getColumn());
            }
        }

        // learning rate for mean of population size variational posterior, default 10000
        ParamExtractor psmeanlrParam = new ParamExtractor("psmeanlr", this.params, this.errorDetected);
        if (psmeanlrParam.ContainsSwitch) {
            if (psmeanlrParam.PostSwitchParam != null) {
                try {
                    Utils.POP_SIZE_MEAN_LEARNING_RATE = Double.parseDouble(psmeanlrParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized learning rate for mean of population size posterior " + psmeanlrParam.PostSwitchValue,
                            psmeanlrParam.PostSwitchParam.getLine(), psmeanlrParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -psmeanlr.",
                        psmeanlrParam.SwitchParam.getLine(), psmeanlrParam.SwitchParam.getColumn());
            }
        }

        // learning rate for standard deviation of node height variational posterior, default 500
        ParamExtractor nhsigmalrParam = new ParamExtractor("nhsigmalr", this.params, this.errorDetected);
        if (nhsigmalrParam.ContainsSwitch) {
            if (nhsigmalrParam.PostSwitchParam != null) {
                try {
                    Utils.NODE_HEIGHT_STDDEV_LEARNING_RATE = Double.parseDouble(nhsigmalrParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized learning rate for standard deviation of node height posterior " + nhsigmalrParam.PostSwitchValue,
                            nhsigmalrParam.PostSwitchParam.getLine(), nhsigmalrParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -nhsigmalr.",
                        nhsigmalrParam.SwitchParam.getLine(), nhsigmalrParam.SwitchParam.getColumn());
            }
        }

        // learning rate for standard deviation of population size variational posterior, default 500
        ParamExtractor pssigmalrParam = new ParamExtractor("pssigmalr", this.params, this.errorDetected);
        if (pssigmalrParam.ContainsSwitch) {
            if (pssigmalrParam.PostSwitchParam != null) {
                try {
                    Utils.POP_SIZE_STDDEV_LEARNING_RATE = Double.parseDouble(pssigmalrParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized learning rate for standard deviation of population size posterior " + pssigmalrParam.PostSwitchValue,
                            pssigmalrParam.PostSwitchParam.getLine(), pssigmalrParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -pssigmalr.",
                        pssigmalrParam.SwitchParam.getLine(), pssigmalrParam.SwitchParam.getColumn());
            }
        }

        // minimum value of node height posterior standard deviation, default 10000
        ParamExtractor nhsigmaminParam = new ParamExtractor("nhsigmamin", this.params, this.errorDetected);
        if (nhsigmaminParam.ContainsSwitch) {
            if (nhsigmaminParam.PostSwitchParam != null) {
                try {
                    Utils.NODE_HEIGHT_MIN_STDDEV = Double.parseDouble(nhsigmaminParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized minimum value of node height posterior standard deviation " + nhsigmaminParam.PostSwitchValue,
                            nhsigmaminParam.PostSwitchParam.getLine(), nhsigmaminParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -nhsigmamin.",
                        nhsigmaminParam.SwitchParam.getLine(), nhsigmaminParam.SwitchParam.getColumn());
            }
        }

        // minimum value of population size posterior standard deviation, default 3000
        ParamExtractor pssigmaminParam = new ParamExtractor("pssigmamin", this.params, this.errorDetected);
        if (pssigmaminParam.ContainsSwitch) {
            if (pssigmaminParam.PostSwitchParam != null) {
                try {
                    Utils.POP_SIZE_MIN_STDDEV = Double.parseDouble(pssigmaminParam.PostSwitchValue);
                } catch (NumberFormatException e) {
                    errorDetected.execute("Unrecognized minimum value of population size posterior standard deviation " + pssigmaminParam.PostSwitchValue,
                            pssigmaminParam.PostSwitchParam.getLine(), pssigmaminParam.PostSwitchParam.getColumn());
                }
            } else {
                errorDetected.execute("Expected value after switch -pssigmamin.",
                        pssigmaminParam.SwitchParam.getLine(), pssigmaminParam.SwitchParam.getColumn());
            }
        }

        noError = noError && checkForUnknownSwitches(
                "st", "mu", "rho", "nhsigma", "pssigma",
                "psp",
                "n0", "r", "nb",
                "ns", "niter", "nhmeanlr", "psmeanlr", "nhsigmalr", "pssigmalr", "nhsigmamin", "pssigmamin"
        );

        checkAndSetOutFile(
                stParam, muParam, rhoParam, nhsigmaParam, pssigmaParam,
                pspParam,
                n0Param, rParam, nbParam,
                nsParam, niterParam, nhmeanlrParam, psmeanlrParam, nhsigmalrParam, pssigmalrParam, nhsigmaminParam, pssigmaminParam
        );

        return noError;
    }

    @Override
    protected String produceResult() {
        /*
         * initialize inference
         */
        Utils.MUTATION_RATE = _originalMu * 4 * Utils.N0;
        RecombinationRate recombRate = new RecombinationRate(_recombRate);
        ModelTree initModel = new ModelTree(_startTree, recombRate, true);
        Prior prior = new Prior();
        VariationalInference algo = new VariationalInference(initModel, prior);

        // log time used by algo.run()
        long startTime = System.currentTimeMillis();
        algo.run();
        long endTime = System.currentTimeMillis();
        long totalTimeMillis = endTime - startTime;

        System.out.println("");
        System.out.println("===================================== Final results =====================================");
        VariationalModel posterior = algo.getVariationalPosterior();
        System.out.println("********************* Node heights in postorder traversal *********************");
        int nodeIndex = 1;
        for (VariationalVariable var:posterior.getNodeHeightVariableList()) {
            System.out.println("Node " + nodeIndex + "; Mean: " + var.getMean() + ", Standard deviation: " + var.getStandardDeviation());
            nodeIndex += 1;
        }
        System.out.println("********************* Population sizes in postorder traversal *********************");
        int branchIndex = 1;
        for (VariationalVariable var:posterior.getPopSizeVariableList()) {
            System.out.println("Branch " + branchIndex + "; Mean: " + var.getMean() + ", Standard deviation: " + var.getStandardDeviation());
            branchIndex += 1;
        }

        System.out.println("");
        System.out.println("Total execution time: " + totalTimeMillis / 1000.0 + " s");
        System.out.println("Time used to build HMM: " + Utils.buildingTime / 1000.0 + " s");
        System.out.println("Time used in forward algorithm: " + Utils.likelihoodTime / 1000.0 + " s");

        return "";
    }

}
