package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.IOException;
import java.io.StringReader;
import java.nio.charset.Charset;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

public class MSDataParser {

    List<Tree> result = new ArrayList<>();


    public MSDataParser(List<String> data)
    {

        for (String line : data)
        {
            String[] parts = line.split(Pattern.quote("]"));

            String numberPart = parts[0].replaceAll(Pattern.quote("["),"");
            int count = Integer.parseInt(numberPart);

            Tree tree = parseTree(parts[1]);

            for (int i = 0; i < count; i++) {
                result.add(tree);
            }

        }

    }

    Tree parseTree(String treePart)
    {
        NewickReader reader = new NewickReader(new StringReader(treePart));
        try {
            return reader.readTree();
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    public List<Tree> getTrees()
    {
        return Collections.unmodifiableList(result);
    }

    public static void main(String[] args) throws IOException {
        Path p = FileSystems.getDefault().getPath("NexusDirectory", "KnownGeneTrees/SecretSauceGeneTrees");
        List<String> lines = Files.readAllLines(p, Charset.defaultCharset());

        MSDataParser par = new MSDataParser(lines);

        List<Tree> trees = par.getTrees();

        List<Tree> distinctTrees = new ArrayList<>();

        double[] treeProportions = {0,0};


        for (Tree t: trees)
        {
            int nodeX = Integer.parseInt(t.toNewick().substring(1, 2));
            if(nodeX==1)
                treeProportions[0] +=1;
            if(nodeX==3)
                treeProportions[1] += 1;
            if(nodeX!=1 && nodeX!=3)
                System.out.println("PROBLEM");



            /*double branchLength = t.getRoot().getChildren().iterator().next().getParentDistance();
            if(branchLength>7.5 && branchLength<10.0)
                System.out.println("Problem..." + branchLength);
            if(branchLength>7.5){
                treeProportions[1] += 1;
            }
            else
                treeProportions[0] += 1;*/


            /*boolean alreadyIn = distinctTrees.stream().anyMatch(oldTree ->
                            Trees.sameTopology(t,oldTree)
            );

            if (!alreadyIn)
                distinctTrees.add(t);*/
        }
        System.out.println("# of Normal Trees = " + treeProportions[0]);
        System.out.println("# of MIG Trees = " + treeProportions[1]);
        System.out.println("% Normal Trees = " + 100*treeProportions[0]/(treeProportions[0]+treeProportions[1]) + "%");
        System.out.println("% MIG Trees = " + 100*treeProportions[1]/(treeProportions[0]+treeProportions[1]) + "%");

        System.out.println(distinctTrees);
    }


}
