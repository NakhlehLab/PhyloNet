package edu.rice.cs.bioinfo.programs.phylonet.commands;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.Program;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.Configuration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.ConfigurationBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.CreateOptimumHMM;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.gui.PlotOutput;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.FullOutputInformation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.HmmOptimizationFunction;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.JsonUtilities;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.ListUtilities;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.MSDataParser;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.ParamParser;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import static edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils.generateParentalTrees;


@CommandName("HmmCommand")
public class HmmCommand extends CommandBaseFileOutWithDNA {
    Map<String,String> sourceIdentToDNA;

    List<Tree> MSTrees;
    Map<String,List<String>> speciesToAlleles;
    Network<String> net;

    boolean plots;

    Path outputDirectory;

    ConfigurationBuilder configBuilder;

    public HmmCommand(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Map<String, String> dna, Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader)
    {
        super(motivatingCommand, params, sourceIdentToNetwork, dna, errorDetected, rnReader);
        sourceIdentToDNA = getSourceIdentToDNA();
    }


    @Override
    protected int getMinNumParams() {
        return 1;
    }

    @Override
    protected int getMaxNumParams() {
        return 20;
    }

    @Override
    protected boolean checkParamsForCommand() {

        ParamParser parser = new ParamParser(params,errorDetected);

        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        net = transformer.makeNetwork(assertAndGetNetwork(0));

        String speciesNetworkString = NetworkTransformer.toENewick(assertAndGetNetwork(0));
        System.out.println(speciesNetworkString);

        speciesToAlleles = parser.getSpeciesToAlleleMapParameter("allelemap");

        if (speciesToAlleles == null)
            speciesToAlleles = getDefaultAlleleMap();

        String treeFile = parser.getStringParameter("treefile");
        if (treeFile != null)
            MSTrees = getTreesFromMSData(treeFile);

        String outputDirectoryPath = parser.getStringParameter("outputdirectory");
        if (outputDirectoryPath != null)
            outputDirectory = Paths.get(outputDirectoryPath);

        plots = !parser.getBooleanParameter("noplots");

        if (!allAllelesInDataSection())
            throw new RuntimeException("Not all of the alleles are in the data section.");

        if (!verifyAllSpeciesInTaxaMap())
        {
            System.out.println("The speciesToAlleles map keys are: " + speciesToAlleles.keySet().toString());
            for(String key: speciesToAlleles.keySet())
            {
                System.out.println("The speciesToAlleles map values are: " + speciesToAlleles.get(key));
            }

            System.out.println("The sourceIdentToDNA map keys are: " + sourceIdentToDNA.keySet().toString());
            throw new RuntimeException("Not all species are in the allele map");
        }

        configBuilder = createConfigurationBuilder(parser);

        parser.checkNoInvalidParameters();

        this.checkAndSetOutFile(parser.getSwitches());

        return true;
    }

    private boolean allAllelesInDataSection()
    {
        for (String species : speciesToAlleles.keySet())
        {
            for (String allele : speciesToAlleles.get(species))
            {
                if (!sourceIdentToDNA.containsKey(allele))
                    return false;
            }
        }
        return true;
    }

    private ConfigurationBuilder createConfigurationBuilder(ParamParser parser)
    {
        ConfigurationBuilder configBuilder;

        configBuilder = ConfigurationBuilder.getNormal();

        Integer threads = parser.getIntegerParameter("threads");
        if (threads != null)
            configBuilder.withThreads(threads);

        Integer iterations = parser.getIntegerParameter("iterations");
        if (iterations != null)
            configBuilder.withIterations(iterations);

        Integer numberOfFolds = parser.getIntegerParameter("ppatraining");
        if (numberOfFolds != null)
            configBuilder.withPPATraining(numberOfFolds);

        Integer numberOfRuns = parser.getIntegerParameter("numberofruns");
        if (numberOfRuns != null)
            configBuilder.withNumberOfRuns(numberOfRuns);

        if (parser.getBooleanParameter("gtr"))
            configBuilder.withGTR();

        if (numberOfFolds != null && numberOfFolds == 1)
        {
            throw new RuntimeException("The value after the -ppatraining switch in the Nexus file is K: the number of " +
                    "folds in the K-Folds training procedure. K is currently equal to 1, which contradicts the purpose of " +
                    "the K-folds training procedure. Try running a standard run by deleting the -ppatraining switch and " +
                    "following value completely, or by setting the value after the switch to an integer greater than 1.");
        }

        Double threshold = parser.getDoubleParameter("threshold");
        if(threshold != null)
            configBuilder.withThreshold(threshold);

        return configBuilder;
    }

    private Map<String, List<String>> getDefaultAlleleMap()
    {
        Map<String, List<String>> result = new HashMap<>();

        for (String allele : sourceIdentToDNA.keySet())
        {
            result.put(allele,Arrays.asList(allele));
        }
        return result;
    }

    private List<Tree> getTreesFromMSData(String MSDataFileName)
    {
        Path inputFile = Paths.get(Program.inputNexusFileName);

        Path MSDataFile = inputFile.resolveSibling(MSDataFileName);

        List<String> geneTreeData;
        try {
            geneTreeData = Files.readAllLines(MSDataFile, Charset.defaultCharset());

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        MSDataParser parseMyData = new MSDataParser(geneTreeData);
        return parseMyData.getTrees();
    }

    private boolean verifyAllSpeciesInTaxaMap()
    {
        for (NetNode<String> node : net.getLeaves())
            if (!speciesToAlleles.containsKey(node.getName()))
            {
                System.out.println(node.getName());
                return false;
            }

        return true;
    }

    @Override
    protected String produceResult() {
        Configuration finalConfig = configBuilder.build(Program.inputNexusFileName);
        List<STITree<String>> sourceTrees =  generateParentalTrees(net, speciesToAlleles);
        System.out.println("\nThe source trees were: "+sourceTrees);


        if(finalConfig.PPATraining)
        {
            runPPATraining(finalConfig, sourceTrees);
        }
        else
        {
            runNormalTraining(finalConfig, sourceTrees);
        }


        return "";
    }

    private FullOutputInformation runPPATraining(Configuration finalConfig, List<STITree<String>> sourceTrees)
    {
        //Declare the 'best-of' structure, which holds the best trained model for each fold, and its training output data.
        List<FullOutputInformation> bestModelsAndTrainingData = new ArrayList<>();

        //Splits the data into K folds.
        List<NucleotideObservation> allDataObservations = NucleotideObservation.getObservations(sourceIdentToDNA);
        List<List<NucleotideObservation>> ObservationFolds = ListUtilities.partition(allDataObservations, finalConfig.numberOfFolds);

        for(int fold=0; fold<finalConfig.numberOfFolds; fold++)
        {
            Path resolvedDirectory = null;
            if (outputDirectory != null)
                resolvedDirectory = outputDirectory.resolve("fold" + fold);

            //Gets the training and testing sets for this fold, as well as the MS Trees fold.
            List<NucleotideObservation> trainingSetObs = ListUtilities.pieceTogether(ObservationFolds, fold);
            List<NucleotideObservation> testingSetObs = ObservationFolds.get(fold);

            //Gets the best trained HMM on this fold.
            bestModelsAndTrainingData.add(getHmm(finalConfig, trainingSetObs, testingSetObs, resolvedDirectory, sourceTrees , null));
        }

        //Gets the best trained model, and writes its training output to a file.
        FullOutputInformation grossBestModel = getBestModel(bestModelsAndTrainingData);

        if (outputDirectory != null)
        {
            Path bestTrainingDirectory = outputDirectory.resolve("overallbesttraining");
            writeOutput(grossBestModel, bestTrainingDirectory);
        }


        //Run the best trained model on the entire data (not just the training data), write output to a file, plot the output.
        FullOutputInformation outputInformation = reprocessOutputInformation(grossBestModel, finalConfig, sourceTrees, allDataObservations);

        if (outputDirectory != null)
        {
            Path bestTestingDirectory = outputDirectory.resolve("overallbesttesting");
            writeOutput(outputInformation, bestTestingDirectory);
        }


        if (plots) {
            new PlotOutput(outputInformation,finalConfig).show();
        }
        return grossBestModel;
    }

    private FullOutputInformation reprocessOutputInformation(FullOutputInformation grossBestModel, Configuration finalConfig, List<STITree<String>> sourceTrees, List<NucleotideObservation> allDataObservations)
    {
        try (HmmOptimizationFunction copyBestModel = new HmmOptimizationFunction(net, speciesToAlleles, allDataObservations, finalConfig))
        {
            copyBestModel.value(grossBestModel.rawOutputInformation.bestParametersRaw);
            return new FullOutputInformation(copyBestModel,MSTrees,sourceTrees,allDataObservations, finalConfig.threshold);
        }
    }

    private FullOutputInformation runNormalTraining(Configuration finalConfig, List<STITree<String>> sourceTrees)
    {
        List<NucleotideObservation> dnaObservations = NucleotideObservation.getObservations(sourceIdentToDNA);

        FullOutputInformation bestRunData = getHmm(finalConfig,
                dnaObservations, dnaObservations, outputDirectory, sourceTrees, MSTrees);

        if (plots)
            new PlotOutput(bestRunData, finalConfig).show();

        return bestRunData;
    }

    private FullOutputInformation getHmm(Configuration config,
                                         List<NucleotideObservation> trainingObservations, List<NucleotideObservation> testingObservations, Path outputDirectory,
                                         List<STITree<String>> sourceTrees, List<Tree> currentMSTrees)
    {

        FullOutputInformation bestRun = null;
        for (int i = 0; i < config.numberOfRuns;i++)
        {
            CreateOptimumHMM untrainedHMM = getUntrainedHMM(config,trainingObservations);

            try (HmmOptimizationFunction solvedModel = untrainedHMM.solve())
            {
                FullOutputInformation result =  new FullOutputInformation(solvedModel,currentMSTrees,sourceTrees,testingObservations, config.threshold);

                String runDirectory = "run" + i;
                if (outputDirectory != null)
                    writeOutput(result,outputDirectory.resolve(runDirectory));

                if (bestRun == null || result.compareTo(bestRun) >0)
                    bestRun = result;
            }
        }


        if (outputDirectory != null)
        {
            writeConfiguration(config, outputDirectory);
            writeOutput(bestRun, outputDirectory.resolve("bestrun"));
        }

        return bestRun;
    }


    private CreateOptimumHMM getUntrainedHMM(Configuration config, List<NucleotideObservation> trainingObservations)
    {
        config.seed = ConfigurationBuilder.createSeed();
        return new CreateOptimumHMM(net, speciesToAlleles, trainingObservations, config);
    }

    private FullOutputInformation getBestModel(
            List<FullOutputInformation> allRuns)
    {
        return Collections.max(allRuns);
    }

    private void writeOutput(FullOutputInformation outputInformation, Path outputDir)
    {
        try
        {
            Files.createDirectories(outputDir);
            Gson gson = new GsonBuilder().setPrettyPrinting().create();

            Path rawOutputFile = outputDir.resolve("rawOutput.json");
            JsonUtilities.writeJson(rawOutputFile, gson, outputInformation.rawOutputInformation);

            Path outputFile = outputDir.resolve("output.json");
            JsonUtilities.writeJson(outputFile, gson, outputInformation.processedOutputInformation);

        } catch (IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    private void writeConfiguration(Configuration config, Path outputDir)
    {
        Gson gson = new GsonBuilder().setPrettyPrinting().create();
        Path configFile = outputDir.resolve("config.json");
        try {
            Files.createDirectories(outputDir);
            JsonUtilities.writeJson(configFile, gson, config);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}