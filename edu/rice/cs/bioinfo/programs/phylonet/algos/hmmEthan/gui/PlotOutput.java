package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.gui;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.Configuration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.FullOutputInformation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.ProcessedOutputInformation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.RawOutputInformation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.SpeciesComparision;

import java.util.Arrays;

public class PlotOutput {

    RawOutputInformation rawOut;
    ProcessedOutputInformation out;
    Configuration config;

    public PlotOutput(FullOutputInformation outputInformation, Configuration config) {
        this.rawOut = outputInformation.rawOutputInformation;
        this.out = outputInformation.processedOutputInformation;
        this.config = config;
    }


    public void show() {
        //Create the output dialog and plots.
        DialogWithTabs dialogBox = new DialogWithTabs();

        addSpeciesTreeTabs(dialogBox);

        if (rawOut.actualSpeciesTrees != null)
            addActualSpeciesTreesTabs(dialogBox);

        if (config.ALGORITHM != Configuration.AlgorithmType.SNAPP)
            addGeneTreeTabs(dialogBox);

        addOptimizationDetailTabs(dialogBox);

        dialogBox.createAndShowGUI();
    }

    private void addSpeciesTreeTabs(DialogWithTabs dialogBox)
    {
        PlotAreaGraph posteriorDecodingSpeciesTreeProbGraph = new PlotAreaGraph("Posterior Decoding Species Tree Probabilities", rawOut.posteriorProbabilityOfSpeciesTrees);
        dialogBox.addPanel("Posterior Decoding: Spec Trees", posteriorDecodingSpeciesTreeProbGraph);

        PlotBarGraph mostLikelySpeciesTreesPlot = new PlotBarGraph("Most Likely Species Tree by Plurality", "Loci of DNA", "Inferred Species Tree by Locus",
                out.mostLikelySpeciesTrees);
        dialogBox.addPanel("Most Likely: Spec Trees",mostLikelySpeciesTreesPlot);

        if(config.threshold != 0.0)
        {
            PlotBarGraph thresholdedSpeciesTreeGraph = new PlotBarGraph("Most Likely Species Trees: With Threshold", "Loci of DNA",
                    "Inferred Species Tree if Prob Above Threshold", out.thresholdedMostLikelySpeciesTrees);
            dialogBox.addPanel("Most Likely Thresholded: Spec Trees", thresholdedSpeciesTreeGraph);
        }

        PlotBarGraph viterbiSpeciesTreePlot = new PlotBarGraph("Viterbi Sequence: Species Trees", "Loci of DNA", "Inferred Species Tree by Locus",
                rawOut.viterbiSpeciesTrees);
        dialogBox.addPanel("Viterbi: Spec Trees", viterbiSpeciesTreePlot);
    }

    private void addGeneTreeTabs(DialogWithTabs dialogBox)
    {
        PlotAreaGraph posteriorDecodingGeneTreeProbGraph = new PlotAreaGraph("Posterior Decoding Probabilities of Gene Trees", rawOut.posteriorProbabilityOfGeneTrees);
        dialogBox.addPanel("Posterior Decoding: Gene Trees", posteriorDecodingGeneTreeProbGraph);


        PlotBarGraph viterbiGeneTreePlot = new PlotBarGraph("Viterbi Sequence: Gene Trees", "Loci of DNA", "Inferred Gene Tree by Locus",
                rawOut.viterbiGeneTrees);
        dialogBox.addPanel("Viterbi: Gene Trees", viterbiGeneTreePlot);
    }

    private void addOptimizationDetailTabs(DialogWithTabs dialogBox)
    {
        PlotLineGraph probScoreLineGraph = new PlotLineGraph("Likelihood of Data given Model: by Iteration", "Iteration",
                "Likelihood Score of Data given Model", rawOut.probScoreData);
        dialogBox.addPanel("Log(P(Data | Model))", probScoreLineGraph);

        String outputInformation = getOutputInformation();
        TextOutputTab outputInformationTab = new TextOutputTab(outputInformation);
        dialogBox.addPanel("Output Text Info", outputInformationTab);

        String inputInfoString = getInputInformation();
        TextOutputTab inputInformationTab = new TextOutputTab(inputInfoString);
        dialogBox.addPanel("Input Text Info", inputInformationTab);
    }

    private String getInputInformation()
    {
        Gson prettyJson = new GsonBuilder().setPrettyPrinting().create();

        String inputInfoString = "~SETTINGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" +
                                            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+"\n";
        inputInfoString += prettyJson.toJson(config)
                        +  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

        return inputInfoString;
    }

    private String getOutputInformation()
    {
        String textTabString = "";
        //Create a string of the data, for the text output tab.
        textTabString += "Random number seed: " + config.seed + "\n\n";
        textTabString += "All HMMData Parameters: " + rawOut.bestParameters + "\n\n";
        textTabString += "P(Observations|Model) 'Score' = " + rawOut.finalScore + "\n\n";
        textTabString += "P(Observations|Actual) 'Score' = " + rawOut.scoreOfActual + "\n\n";
        textTabString += "Probability of remaining in Species Tree: " + Arrays.toString(rawOut.bestParameters.speciesStayProbs) + "\n";
        textTabString += "Probability of remaining in Gene Tree: " + rawOut.bestParameters.geneStayProb + "\n\n";
        textTabString += "HMM State Definitions: \n";
        textTabString += rawOut.formattedHmmStates + "\n\n";
        textTabString += Arrays.toString(rawOut.formattedHmmStateCounts) + "\n\n";

        textTabString += "\n Species Tree Indices Definitions: Apply to Legends of Plots and Output Text Info \n";

        String formatSpeciesTrees = "";
        for (int i  =0; i < rawOut.formattedSpeciesTrees.length; i++)
            formatSpeciesTrees += String.format("%d: %s\n",i,rawOut.formattedSpeciesTrees[i]);
        textTabString += formatSpeciesTrees+ "\n\n";

        textTabString += "\n\n" + rawOut.formattedNetwork + "\n\n";

        textTabString += new GsonBuilder().setPrettyPrinting().create().toJson(out);
        return textTabString;
    }

    private void addActualSpeciesTreesTabs(DialogWithTabs dialogBox)
    {
        PlotBarGraph actualSpeciesTreePlot = new PlotBarGraph("Reality Species Trees: Constructed from MS Gene Trees", "Loci of DNA", "Actual Species Tree by Locus", rawOut.actualSpeciesTrees);

        dialogBox.addPanel("Reality: Spec Trees", actualSpeciesTreePlot);


        System.out.println("Viterbi species topologies percent error = " + out.percentErrorBetweenActualAndViterbi);

        dialogBox.addPanel("Viterbi & Reality Comparison: Spec Trees", SpeciesComparision.getCumulativeErrorPlot(rawOut.viterbiSpeciesTrees, rawOut.actualSpeciesTrees));
    }
}
