package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.util.CombinatoricsUtils;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 10/18/18
 * Time: 2:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class LizardDataset {
    List<String> sampleNames = new ArrayList<>();
    Multimap<String, String> species2samples = null;
    Map<String, List<String>> species2alleles = new HashMap<>();
    Map<String, String> alleles2species = new HashMap<>();
    Map<String, Map<String, String>> loci2seq = new HashMap<>();
    int nchar = 0;
    boolean sgt = true;
    String outgroup = "Lampropholis_guichenoti";
    String gtoutgroup = "SP07_indexing28_h0";
    String root;

    void input1() {
        species2samples.put("coggeri", "SP02A_indexing4");
        species2samples.put("coggeri", "SP02A_indexing5");

        species2samples.put("timlowi", "SP09_indexing21");
        species2samples.put("timlowi", "SP09_indexing22");

        species2samples.put("amax", "SP03_indexing25");
        species2samples.put("amax", "SP03_indexing26");

        species2samples.put("insularis", "AS01_indexing45");
        species2samples.put("insularis", "AS01_indexing46");

        species2samples.put("johnstonei", "AS01_indexing29");
        species2samples.put("johnstonei", "AS01_indexing30");


        species2samples.put("triacantha", "AS01_indexing34");
        species2samples.put("triacantha", "AS01_indexing39");

        species2samples.put("tetradactyla", "SP07_indexing20");
        species2samples.put("tetradactyla", "SP08_indexing20");

        species2samples.put("rufilatus", "SP04_indexing9");
        species2samples.put("rufilatus", "SP05_indexing16");

        species2samples.put("munda", "SP04_indexing43");
        species2samples.put("munda", "SP07_indexing10");

        species2samples.put("gracilis", "SP05_indexing14");
        species2samples.put("gracilis", "SP05_indexing2");

        species2samples.put("jarnoldae", "SP07_indexing8");
        species2samples.put("jarnoldae", "SP08_indexing8");

        species2samples.put("rubigo", "SP07_indexing17");
        species2samples.put("rubigo", "SP08_indexing17");

        species2samples.put("pectoralis", "SP07_indexing14");
        species2samples.put("pectoralis", "SP08_indexing14");

        species2samples.put("inconnexa", "SP08_indexing7");
        species2samples.put("inconnexa", "SP10_indexing56");

        species2samples.put("dogare", "SP07_indexing5");
        species2samples.put("dogare", "SP08_indexing5");

        species2samples.put("vivax", "SP09_indexing3");
        species2samples.put("vivax", "SP09_indexing4");

        species2samples.put("decora", "SP07_indexing4");
        species2samples.put("decora", "SP08_indexing4");

        species2samples.put("mundivensis", "SP07_indexing11");
        species2samples.put("mundivensis", "SP08_indexing11");
    }

    void input2() {

        species2samples.put("Liburnascincus_mundivensis", "SP07_indexing11");
        species2samples.put("Liburnascincus_coensis", "SP09_indexing17");
        species2samples.put("Liburnascincus_artemis", "SP09_indexing29");
        species2samples.put("Liburnascincus_scirtetis", "SP09_indexing19");
        species2samples.put("Lygisaurus_macfarlani", "SP09_indexing39");
        species2samples.put("Lygisaurus_sesbrauna", "SP09_indexing15");
        species2samples.put("Lygisaurus_parrhasius", "SP07_indexing13");
        species2samples.put("Lygisaurus_foliorum", "SP07_indexing29");
        species2samples.put("Lygisaurus_aeratus", "SP09_indexing5");
        species2samples.put("Carlia_longipes", "SP07_indexing9");
        species2samples.put("Carlia_rhomboidalis", "SP10_indexing20");
        species2samples.put("Carlia_vivax", "SP09_indexing3");
        species2samples.put("Carlia_rubigo", "SP07_indexing17");
        species2samples.put("Carlia_jarnoldae", "SP07_indexing8");
        species2samples.put("Carlia_amax", "SP03_indexing25");
        species2samples.put("Pygmaeascincus_timlowi", "SP09_indexing21");
        species2samples.put("Lampropholis_guichenoti", "SP07_indexing28");
        species2samples.put("Lampropholis_coggeri", "SP02A_indexing4");
    }

    void input11taxa() {

        species2samples.put("Liburnascincus_mundivensis", "SP07_indexing11");
        //species2samples.put("Liburnascincus_coensis", "SP09_indexing17");
        //species2samples.put("Liburnascincus_artemis", "SP09_indexing29");
        //species2samples.put("Liburnascincus_scirtetis", "SP09_indexing19");
        species2samples.put("Lygisaurus_macfarlani", "SP09_indexing39");
        species2samples.put("Lygisaurus_sesbrauna", "SP09_indexing15");
        //species2samples.put("Lygisaurus_parrhasius", "SP07_indexing13");
        species2samples.put("Lygisaurus_foliorum", "SP07_indexing29");
        //species2samples.put("Lygisaurus_aeratus", "SP09_indexing5");
        species2samples.put("Carlia_longipes", "SP07_indexing9");
        species2samples.put("Carlia_rhomboidalis", "SP10_indexing20");
        species2samples.put("Carlia_vivax", "SP09_indexing3");
        //species2samples.put("Carlia_rubigo", "SP07_indexing17");
        //species2samples.put("Carlia_jarnoldae", "SP07_indexing8");
        species2samples.put("Carlia_amax", "SP03_indexing25");
        species2samples.put("Pygmaeascincus_timlowi", "SP09_indexing21");
        species2samples.put("Lampropholis_guichenoti", "SP07_indexing28");
        species2samples.put("Lampropholis_coggeri", "SP02A_indexing4");
    }

    void input4taxa() {

        species2samples.put("Carlia_rhomboidalis", "SP10_indexing20");
        species2samples.put("Carlia_vivax", "SP09_indexing3");
        species2samples.put("Carlia_amax", "SP03_indexing25");
        species2samples.put("Lampropholis_guichenoti", "SP07_indexing28");
    }

    void input3taxa() {

        species2samples.put("Carlia_longipes", "SP07_indexing9");
        species2samples.put("Carlia_vivax", "SP09_indexing3");
        species2samples.put("Carlia_amax", "SP03_indexing25");
    }

    LizardDataset(String root) {
        this.root = root;

        species2samples = HashMultimap.create();
        input11taxa();


        for(String species : species2samples.keySet()) {
            species2alleles.put(species, new ArrayList<>());

            for(String sample : species2samples.get(species)) {
                species2alleles.get(species).add(sample + "_h0");
                species2alleles.get(species).add(sample + "_h1");
            }
        }

        for(String species : species2alleles.keySet()) {
            for(String allele : species2alleles.get(species)) {
                alleles2species.put(allele, species);
            }
        }


    }

    static void writeFastaFile(Map<String, String> seq, String filename) {
        try {
            PrintWriter out = new PrintWriter(filename);

            for(String taxon : seq.keySet()) {
                out.println(">" + taxon);
                String cur = seq.get(taxon);
                for(int begin = 0, end = Math.min(80, cur.length()) ; begin < cur.length() ; ) {
                    out.println(cur.substring(begin, end));
                    begin = end;
                    end = Math.min(end + 80, cur.length());
                }
                out.println();
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static Map<String, String> readPhyFile(String filename) {
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String s;
            int index = 0;
            s = in.readLine();
            s = s.trim();
            String ss[] = s.split("\\s+");
            int n = Integer.parseInt(ss[0]); // rows
            int m = Integer.parseInt(ss[1]); // colums
            List<String> taxonNames = new ArrayList<>();
            Map<String, String> sequences = new HashMap<>();

            while((s = in.readLine()) != null) {
                s = s.trim();
                if(s.length() < 2) continue;
                if(taxonNames.size() < n) {
                    ss = s.split("\\s+");
                    taxonNames.add(ss[0]);
                    sequences.put(ss[0], "");
                } else {
                    ss[0] = taxonNames.get(index);
                    ss[1] = s;
                }
                sequences.put(ss[0], sequences.get(ss[0]) + ss[1]);

                index++;
                if(index >= n) index = 0;
            }

            for(String taxon : sequences.keySet()) {
                if(sequences.get(taxon).length() != m) {
                    return null;
                }
            }

            return sequences;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    static List<String> readLociList(String filename) {
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String s;
            int index = 0;
            List<String> result = new ArrayList<>();

            while((s = in.readLine()) != null) {
                s = s.trim();

                result.add(s);
                index++;
            }

            System.out.println(result.size());

            return result;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    static Map<String, String> readIQTree(String resultFolder) {
        File path = new File(resultFolder);

        List<String> filenames = new ArrayList<>();
        Map<String, String> filename2locusname = new HashMap<>();
        Map<String, String> locus2treestring = new TreeMap<>();

        File[] files = path.listFiles();

        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){
                if(files[i].getName().endsWith(".treefile")) {
                    String ss[] = files[i].getName().split("_");
                    filename2locusname.put(files[i].toString(), ss[1] + "_" + ss[2]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);
        System.out.println("IQTREE read: " + filenames.size());
        for(String filename : filenames) {
            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                s = in.readLine();
                Tree tree = Trees.readTree(s);
                locus2treestring.put(filename2locusname.get(filename), s);

            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        return locus2treestring;
    }

    void writeToNEXUS() {

        int total = species2alleles.size();
        if (outgroup != null) total--;

        total = (int) CombinatoricsUtils.binomialCoefficient(total, 3);
        System.out.println("Subproblems: " + total);

        List<Integer> indices = new ArrayList<>();
        for (int nnn = 0; nnn < total; nnn++) {
            indices.add(nnn);
        }
        Collections.shuffle(indices, new Random(12345678L));
        Map<String, String> locus2tree = null;
        if(sgt) {
            locus2tree = readIQTree(root + "/data/iqtree/");
        }

        for(int nnn = 0 ; nnn < total ; nnn++) {
            String nexusPath = root + "/run_" + nnn + ".nex";
            System.out.println(nexusPath);

            try {
                PrintWriter out = new PrintWriter(nexusPath);

                out.print(String.format("#NEXUS \nBegin data;\nDimensions ntax=%d nchar=%d;\nFormat datatype=dna symbols=\"012\" missing=? gap=-;\nMatrix\n\n", species2alleles.size(), nchar));

                for (String locus : loci2seq.keySet()) {
                    int n = loci2seq.get(locus).values().iterator().next().length();
                    out.println(String.format("[%s, %d, ...]", locus, n));
                    for (String allele : loci2seq.get(locus).keySet()) {
                        out.println(String.format("%s %s", allele, loci2seq.get(locus).get(allele)));
                    }
                }

                out.println();
                out.println(";End;\n");

                if(sgt) {
                    out.println("BEGIN TREES;");
                    for(String locus : locus2tree.keySet()) {
                        out.println("Tree " + locus + " = " + locus2tree.get(locus));
                    }
                    out.println("End;");
                    out.println();
                }

                out.println("BEGIN PHYLONET; \n");
                out.println("SN_SEQ -cl 6000000 -bl 3000000 -sf 5000 -ee -sd 123456 -subsetloci 150 -fixgttopo -gtburnin -pre 10 -varyps -pl 2 ");
                if(outgroup != null) {
                    out.println("-outgroup \"" + outgroup +"\"");
                }

                if(gtoutgroup != null) {
                    out.println("-gtoutgroup \"" + gtoutgroup +"\"");
                }

                if(sgt) {
                    out.print("-sgt (");
                    boolean first = true;
                    for(String locus : locus2tree.keySet()) {
                        if(first) first = false;
                        else out.print(",");
                        out.print("" + locus);
                    }
                    out.println(") ");
                }

                out.print("-tm <");
                int index = 0;
                for (String species : species2alleles.keySet()) {
                    if (index > 0) out.print(";");
                    out.print(species + ":");
                    for (int i = 0; i < species2alleles.get(species).size(); i++) {
                        if (i > 0) out.print(",");
                        out.print(species2alleles.get(species).get(i));
                    }
                    index++;
                }

                out.println("> ");
                out.println(" -nnn " + nnn + " ;");
                out.println("END;");

                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

             //break;
        }
    }

    void writeFullToNEXUS(String nexusPath, String lociListPath) {

        List<String> lociList = readLociList(lociListPath);

        try {
            PrintWriter out = new PrintWriter(nexusPath);

            out.print(String.format("#NEXUS \nBegin data;\nDimensions ntax=%d nchar=%d;\nFormat datatype=dna symbols=\"012\" missing=? gap=-;\nMatrix\n\n", species2alleles.size(), nchar));

            for (String locus : loci2seq.keySet()) {
                if(!lociList.contains(locus)) continue;

                int n = loci2seq.get(locus).values().iterator().next().length();
                out.println(String.format("[%s, %d, ...]", locus, n));
                for (String allele : loci2seq.get(locus).keySet()) {
                    out.println(String.format("%s %s", allele, loci2seq.get(locus).get(allele)));
                }
            }

            out.println();
            out.println(";End;\nBEGIN PHYLONET; \n");
            out.println("MCMC_SEQ -cl 10000000 -bl 3000000 -sf 5000 -ee -sd 123456 -pl 4 ");
            out.print("-tm <");
            int index = 0;
            for (String species : species2alleles.keySet()) {
                if (index > 0) out.print(";");
                out.print(species + ":");
                for (int i = 0; i < species2alleles.get(species).size(); i++) {
                    if (i > 0) out.print(",");
                    out.print(species2alleles.get(species).get(i));
                }
                index++;
            }

            out.println("> ;");
            out.println("END;");

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void read() {
        List<String> loci = readLociList(root + "/data/lists/loci_complete_informative_data_filtered.txt");

        for(String locus : loci) {
            String locusPath = root + "/data/phy/" + locus + ".NT.TNR.trim.gt.phy";
            Map<String, String> onelocus = readPhyFile(locusPath);
            Map<String, String> selectedAlleles = new HashMap<>();
            for(String allele : alleles2species.keySet()) {
                selectedAlleles.put(allele, onelocus.get(allele));
            }
            loci2seq.put(locus, selectedAlleles);
        }
    }

    static void prepare(String[] args) {
        String path = args.length > 0 ? args[0] : "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/lizard";
        LizardDataset lizardDataset = new LizardDataset(path);

        lizardDataset.read();

        lizardDataset.writeToNEXUS();

    }

    static void prepare_astral(String[] args) {
        LizardDataset lizardDataset = new LizardDataset(args[0]);
        lizardDataset.read();

        try {


            PrintWriter out = new PrintWriter(args[0] + "/data/astral/astral_input.txt");

            for(String locus : lizardDataset.loci2seq.keySet()) {
                BufferedReader in = new BufferedReader(new FileReader(args[0] + "/data/iqtree/" + locus + ".NT.TNR.trim.gt.phy.treefile"));
                String s = in.readLine();
                out.println(s);
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    static void prepare_iqtree(String[] args) {
        LizardDataset lizardDataset = new LizardDataset(args[0]);
        lizardDataset.read();

        for(String locus : lizardDataset.loci2seq.keySet()) {

            writeFastaFile(lizardDataset.loci2seq.get(locus), args[0] + "/data/iqtree_4taxa/iqtree_" + locus + "_input.fasta");
        }



    }

    static void check() {
        String resultFolder = "/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/results/run70";
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
        options.outgroup = "Lampropholis_guichenoti"; //"Lampropholis_coggeri";//"Pygmaeascincus_timlowi"; //Lampropholis_guichenoti
        options.eps = 0.01;
        //options.trueNetwork = Networks.readNetwork("((((Lampropholis_coggeri:0.005421004068874944)I3#H1:0.009292050652544074,((Lygisaurus_macfarlani:0.0032689514794862575)I7#H2:0.007939200985418114,((((Carlia_vivax:0.003990805427495498)I11#H3:0.0036997156781629396,(Carlia_amax:0.002415459661965271)I12#H4:0.005275061443693167)I9:0.0021545692750479044,((I12#H4:0.005595162676113067,((Carlia_longipes:0.005875135856075196,Carlia_rhomboidalis:0.005875135856075196)I17:3.257814894150219E-4,I11#H3:0.0022101119179947218)I15:0.0018097049925881181)I13:0.0014377306412839724,(Lygisaurus_foliorum:0.007092727836841509,(Lygisaurus_sesbrauna:0.005873230850686941,I7#H2:0.0026042793712006836)I16:0.0012194969861545673)I14:0.0023556251425208012)I10:3.96737401344032E-4)I8:5.545968113897015E-4,Liburnascincus_mundivensis:0.010399687192096042)I6:8.084652728083285E-4)I5:0.003504902256514646)I4:7.926837682867198E-4,Pygmaeascincus_timlowi:0.015505738489705736)I2:0.008717091433705469,(I3#H1:0.007819181189534181,Lampropholis_guichenoti:0.013240185258409125)I1:0.01098264466500208)I0;");

        SuperNetwork3.printDetails_ = true;
        //SNSummary summary = Pipeline.stage2(filenames, 6000000, 3000000, 5000, options);
        SNSummary summary = Pipeline.stage2_1(filenames, 1000000, 500000, 5000, options);

        Network inferred = summary.inferredNetwork;
        System.out.println(inferred);
    }

    static void printAstralMapping() {
        LizardDataset lizardDataset = new LizardDataset(".");
        System.out.println();
        boolean first0 = true;
        for(String species : lizardDataset.species2alleles.keySet()) {

            System.out.print(species);
            System.out.print(":");
            boolean first1 = true;
            for(String allele : lizardDataset.species2alleles.get(species)) {
                if(first1) first1 = false;
                else System.out.print(",");
                System.out.print(allele);
            }
            System.out.println();

        }
    }

    static void prepareMCMCSEQ() {
        LizardDataset lizardDataset = new LizardDataset("/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/lizard");
        lizardDataset.read();
        lizardDataset.writeFullToNEXUS("/Users/zhujiafan/Documents/BioinfoData/Debug/mcmcseq.nex", "/Users/zhujiafan/Documents/BioinfoData/Debug/locilist.txt");
    }

    static void printMapping() {
        LizardDataset lizardDataset = new LizardDataset(".");
        System.out.println();
        System.out.print("<");
        boolean first0 = true;
        for(String species : lizardDataset.species2alleles.keySet()) {
            if(first0) first0 = false;
            else System.out.print(";");

            System.out.print(species);
            System.out.print(":");
            boolean first1 = true;
            for(String allele : lizardDataset.species2alleles.get(species)) {
                if(first1) first1 = false;
                else System.out.print(",");
                System.out.print(allele);
            }

        }
        System.out.print(">");
    }

    public static void main(String[] args) {
        //prepareMCMCSEQ();
        //printMapping();
        //printAstralMapping();
        //prepare_astral(new String[]{"/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/lizard/"});
        //prepare_iqtree(new String[]{"/Users/zhujiafan/Documents/BioinfoData/SuperNetwork/lizard/"});
        //prepare(args);
        check();


    }
}
