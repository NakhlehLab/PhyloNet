package edu.rice.cs.bioinfo.programs.vaal2sequencings;

import edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.ir.LineDifferenceRecord;
import edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.parsers.VAALOutStringSplitParser;
import edu.rice.cs.bioinfo.library.language.vaal.out.reading.OutReadResult;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 3/21/13
 * Time: 3:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program
{
    public static void main(String[] args) throws Exception
    {
        File inFile = new File(args[0]);
        Map<File,String> outFileToSourceId = readInFile(inFile);

        VAALOutStringSplitParser outFileParser = new VAALOutStringSplitParser();

        Map<NucleotideLocation, String> referenceGenome = new HashMap<NucleotideLocation, String>();
        Map<File, Map<NucleotideLocation,String>> outFileToDifferences = new LinkedHashMap<File, Map<NucleotideLocation,String>>();
        for(File outFile : outFileToSourceId.keySet())
        {
            String sourceId = outFileToSourceId.get(outFile);
            OutReadResult<List<LineDifferenceRecord>> outContents;
            try
            {
                outContents = outFileParser.parse(new FileInputStream(outFile));
            }
            catch (FileNotFoundException e)
            {
                String path = inFile.getParentFile().getAbsolutePath() + File.separator + outFile.getName();
                System.out.println(path);
                outContents = outFileParser.parse(new FileInputStream(path));
            }

            outFileToDifferences.put(outFile, new HashMap<NucleotideLocation, String>());

            for(LineDifferenceRecord difference : outContents.getDifferences())
            {
                NucleotideLocation differencePosition = new NucleotideLocation(difference.Contig, difference.Position);
                if(!referenceGenome.containsKey(differencePosition))
                {
                    referenceGenome.put(differencePosition, difference.Ref);
                }

                outFileToDifferences.get(outFile).put(differencePosition, difference.Sample);
            }
        }

        List<NucleotideLocation> globalDifferenceLocations = new LinkedList<NucleotideLocation>(referenceGenome.keySet());
        Collections.sort(globalDifferenceLocations);
        System.out.println(globalDifferenceLocations);

        StringBuffer sequencingsText = new StringBuffer("<Sequencings>\n");
        for(File outFile : outFileToDifferences.keySet())
        {
            Map<NucleotideLocation,String> positionToDifference = outFileToDifferences.get(outFile);

            StringBuffer sequenceAttributeValue = new StringBuffer();
            for(NucleotideLocation differenceLocation : globalDifferenceLocations)
            {
                sequenceAttributeValue.append(positionToDifference.containsKey(differenceLocation) ?
                                                positionToDifference.get(differenceLocation) :
                                                referenceGenome.get(differenceLocation));
            }

            sequencingsText.append("\t<Sequencing sourceId=\"" + outFileToSourceId.get(outFile) +"\" ");
            sequencingsText.append("sequence=\"" + sequenceAttributeValue +"\"/>\n");
        }
        sequencingsText.append("</Sequencings>");

        System.out.print(sequencingsText);


    }

    private static Map<File, String> readInFile(File inFile) throws ParserConfigurationException, IOException, SAXException{

        HashMap<File,String> outFileToSourceId = new LinkedHashMap<File, String>();

        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

        // use the factory to take an instance of the document builder
        DocumentBuilder db = dbf.newDocumentBuilder();
        // parse using the builder to get the DOM mapping of the
        // XML file
        Document dom = db.parse(inFile);


        Element doc = dom.getDocumentElement();
        if(doc.getTagName().toLowerCase().equals("vaal2sequencings"))
        {
            NodeList isoloateElements = dom.getElementsByTagName("Isolate");
            for(int i = 0; i<isoloateElements.getLength(); i++)
            {
                Node element = isoloateElements.item(i);
                String outFileText = element.getAttributes().getNamedItem("outFile").getNodeValue();
                String sourceIdText = element.getAttributes().getNamedItem("sourceId").getNodeValue();

                File outFile = new File(outFileText);

                outFileToSourceId.put(outFile, sourceIdText);
            }

            return outFileToSourceId;
        }

        throw new IllegalArgumentException();
    }
}
