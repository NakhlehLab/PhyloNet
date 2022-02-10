package edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.experiments;

import edu.rice.cs.bioinfo.programs.phylonet.algos.dissimilarity.WeightedAveragePathDistance;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class NakhlehMCMCSEQ {
    public static void main(String[] args) {
        if (args.length != 2) {
            System.err.println("Command-line arguments:\n" +
                    "\t[.out file] [results directory]");
            System.exit(-1);
        }

        List<Network<BniNetwork>> networks = new ArrayList<>();


        try {
            Scanner fileReader = new Scanner(new File(args[0]));

            while (fileReader.hasNextLine()) {
                String line = fileReader.nextLine();

                if (line.startsWith("[")) {
                    String richNewick = line.substring(line.indexOf("]") + 1);
                    networks.add(Networks.readNetwork(richNewick));
                }
            }
            fileReader.close();

        } catch (FileNotFoundException e) {
            System.err.println("Could not find .out file to read MCMC_SEQ posterior samples from.");
            System.exit(-1);
        }

        File saveDir = new File(args[1]);

        if (!saveDir.exists() || !saveDir.isDirectory()) {
            System.err.println("The results directory to write experiment results to does not exist.");
            System.exit(-1);
        }

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(Paths.get(args[1], "mcmc_seq.txt").toFile()));

            writer.write("MCMCSEQ experiment\n");
            writer.write("X,Y,Z\n");

            for (int i = 0; i < networks.size(); i++) {
                for (int j = 0; j <= i; j++) {
                    double dissimilarity = Networks.computeDistanceBetweenTwoNetworks(networks.get(i), networks.get(j));
                    writer.write(i + "," + j + "," + dissimilarity + "\n");
                }
                System.out.println((i + 1) + "/" + networks.size());
            }

            writer.close();

        } catch (IOException e) {
            System.err.println("Could not write experiment results.");
            e.printStackTrace();
            System.exit(-1);
        }

    }
}
