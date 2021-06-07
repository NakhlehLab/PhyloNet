package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util;

import com.google.gson.Gson;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;

public class JsonUtilities
{
    public static void writeJson(Path fileName, Gson g, Object obj)
    {
        try (BufferedWriter writer = Files.newBufferedWriter(fileName, Charset.defaultCharset()))
        {
            g.toJson(obj,writer);
        } catch (IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    public static <T> T readJson(Path fileName, Gson g, Class<T> cl)
    {
        try(BufferedReader reader = Files.newBufferedReader(fileName,Charset.defaultCharset()))
        {
            return g.fromJson(reader,cl);

        } catch (IOException e)
        {
            throw new RuntimeException(e);
        }

    }
}
